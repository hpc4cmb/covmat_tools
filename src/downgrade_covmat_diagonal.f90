program degrade_covmat

  ! O'Dwyer 22-04-2010

  USE healpix_types
  USE fitstools
  USE head_fits
  USE pix_tools
  USE extension
  USE udgrade_nr
  USE alm_tools
  USE covmat_util, only : tic, toc

  implicit none

  CHARACTER(len=filenamelen) :: argument, matrix_filein, matrix_fileout
  INTEGER(i4b) :: nsideout, iargument
  real(dp) :: rcondlim
  LOGICAL :: there, invert=.false.

  integer, external :: omp_get_num_procs, omp_get_max_threads, &
       omp_get_thread_num
  integer :: nprocs, nthreads_max, nthreads, id_thread, ierr

  nprocs = omp_get_num_procs()
  nthreads_max = omp_get_max_threads()
  nthreads = nthreads_max
  id_thread = omp_get_thread_num()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in command line arguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (id_thread == 0) &
       write (*,'("OMP: ",i0," procs per node, ",i0, &
       & " threads per task.")') nprocs, nthreads

    IF (nArguments() < 2)  THEN
       write (*,'(/,a,/)') '  Usage: downgrade_covmat_diagonal ' // &
            '<covmat_in> <nside_out> [-o <covmat_out>] [--invert <rcond_threshold>]'
       stop
    ELSE
       CALL getArgument(1, argument)
       matrix_filein = TRIM(argument)
       INQUIRE(file=TRIM(matrix_filein), exist=there)
       IF (.NOT. there) STOP 'covmat_file does not exist!'
       WRITE (*,'(a,i0,a)') &
            'covmat in is == ' // trim(matrix_filein)

       CALL getArgument(2, argument)
       READ (argument, *, iostat=ierr) nSideout
       if (ierr /= 0) stop 'Failed to parse nsideout'
       WRITE (*,'(a,i4)') 'Nsideout   == ', nSideout
       
       iArgument = 3
       matrix_fileout = '!downgraded_' // trim(matrix_filein)
       do
          if (iargument > narguments()) exit
          CALL getArgument(iargument, argument)
          if (index(argument, '-o') > 0) then
             iargument = iargument + 1
             CALL getArgument(iargument, argument)
             matrix_fileout = TRIM(ADJUSTL(argument))
          else if (index(argument, '-inv') > 0) then
             invert = .true.
             iargument = iargument + 1
             CALL getArgument(iargument, argument)
             READ (argument, *, iostat=ierr) rcondlim
             if (ierr /= 0) stop 'Failed to parse rcondlim'
             write (*,*) ' Inverting the covariance matrices and rejecting rcond below ',rcondlim
          else
             write (*,'("ERROR: unrecongnized command line argument: ")') &
                  trim(adjustl(argument))
             stop
          end if
          iargument = iargument + 1
       end do

       WRITE (*,*) ' Writing to ' // TRIM(matrix_fileout)
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL covmat_degrade(matrix_filein, nsideout, matrix_fileout, invert, rcondlim)



CONTAINS



  SUBROUTINE covmat_degrade(matrix_filein, nside_out, matrix_fileout, invert_matrix, rcondlim)

    USE healpix_types
    USE fitstools, only : input_tod, write_bintabh
    USE pix_tools
    USE extension

    IMPLICIT NONE

    character(len=filenamelen) :: matrix_filein, matrix_fileout
    logical,intent(in) :: invert_matrix
    real(dp), intent(in) :: rcondlim
    CHARACTER(len=80) :: header(1000)
    integer(i4b) :: nside_in, nside_out, nmaps, nlheader
    INTEGER(i4b) :: ncols, ierr, ordering
    INTEGER(i8b) :: npix_tot, npix_in, npix_out, ipix_in, ipix_out, npix2one
    REAL(dp), pointer :: cc_one(:, :), cc_sum(:, :)
    REAL(sp), pointer :: cc_out(:, :), cc(:, :)
    ! For inverting the 3x3 matrices using lapack
    real(dp), allocatable :: eigenvals(:), eigenvecs(:,:)
    integer(i8b), parameter :: workspacelen=100
    real(dp), allocatable :: workspace(:)
    integer(i4b) :: info, ieigen

    call tic

    npix_tot = getsize_fits(matrix_filein, nmaps=ncols, ordering=ordering, &
         nside=nside_in)    

    nlheader = 1000
    header = ''
    nmaps = 3
    IF (ncols == 1) nmaps = 1

    npix_in = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)
    if (npix_in < 0 .or. npix_out < 0) then
       write (*,'(" Illegal nside: ",i0,", ",i0)') nside_in, nside_out
       return
    end if
    npix2one = (nside_in / nside_out)**2

    ALLOCATE(cc(0:npix_in-1, ncols), cc_out(0:npix_out-1, ncols), stat=ierr)
    IF (ierr /= 0) STOP 'No room for cc'

    CALL input_tod(matrix_filein, cc, int(npix_in, i8b), ncols, header=header)

    call toc('Read matrix')

    if (ordering==1) then
       call tic
       write (*,*) 'Reordering the matrix into NESTED'
       !call convert_inplace(ring2nest, cc)
       call convert_ring2nest(nside_in, cc)
       call toc('Reorder')
    end if
    
    call add_card(header, 'NSIDE', nside_out)
    call add_card(header, 'LASTPIX', npix_out-1)
    call add_card(header, 'ORDERING', 'NESTED')

    call tic

    !$OMP PARALLEL DEFAULT(PRIVATE) &
    !$OMP   SHARED(nmaps, npix_out, npix2one, cc, cc_out, rcondlim, invert_matrix)
    ALLOCATE(cc_one(nmaps,nmaps), cc_sum(nmaps,nmaps), &
         eigenvecs(nmaps,nmaps), eigenvals(nmaps), workspace(workspacelen))
    !$OMP DO SCHEDULE(STATIC, 4)
    DO ipix_out = 0, npix_out-1
       ! Sum the Nobs matrices of the subpixels, map is NESTED

       cc_sum = 0.0_dp
       DO ipix_in = ipix_out*npix2one, (ipix_out+1)*npix2one-1
          if (cc(ipix_in, 1) == HPX_DBADVAL) cycle
          cc_one(1,  1) = cc(ipix_in, 1)
          IF (nmaps == 3) THEN
             cc_one(1,2:3) = cc(ipix_in, 2:3)
             cc_one(2,  1) = cc(ipix_in, 2)
             cc_one(2,2:3) = cc(ipix_in, 4:5)
             cc_one(3,  1) = cc(ipix_in, 3)
             cc_one(3,  2) = cc(ipix_in, 5)
             cc_one(3,  3) = cc(ipix_in, 6)
          END IF
          cc_sum = cc_sum + cc_one
       END DO

       if (invert_matrix .and. cc_sum(1, 1) /= 0) then
          IF (nmaps == 1) THEN
             cc_sum(1, 1) = 1. / cc_sum(1, 1)
          ELSE
             ! Invert the summed Nobs matrix

             ! compute the eigenvalue decomposition
             ! http://www.netlib.org/lapack/double/dsyev.f
             eigenvecs = cc_sum
             call dsyev('V', 'U', nmaps, eigenvecs, nmaps, eigenvals, &
                  workspace, workspacelen, info)
             if (info == 0 .and. eigenvals(3)/eigenvals(1) > rcondlim) then
                do ieigen = 1,nmaps
                   eigenvecs(:,ieigen) = &
                        eigenvecs(:,ieigen) / sqrt(eigenvals(ieigen))
                end do
                cc_sum = matmul(eigenvecs, transpose(eigenvecs))
             else
                if (info < 0) write (*,'("dsyev failed. info == ",i0)') info
                cc_sum = 0
             end if
          END IF
       end if

       cc_out(ipix_out, 1) = real(cc_sum(1, 1), sp)
       if (nmaps == 3) then
          cc_out(ipix_out, 2) = real(cc_sum(1, 2), sp)
          cc_out(ipix_out, 3) = real(cc_sum(1, 3), sp)
          cc_out(ipix_out, 4) = real(cc_sum(2, 2), sp)
          cc_out(ipix_out, 5) = real(cc_sum(2, 3), sp)
          cc_out(ipix_out, 6) = real(cc_sum(3, 3), sp)
       end if
    END DO
    !$OMP END DO
    deallocate(cc_one, cc_sum, eigenvecs, eigenvals, workspace)
    !$OMP END PARALLEL

    call toc('Downgrade')

    call tic

    ! write out the (degraded) noise matrix here
    CALL write_bintabh(cc_out, int(npix_out, i8b), ncols, header, nlheader, &
         matrix_fileout, repeat=int(npix_out, i4b))

    call toc('Write matrix')

  END SUBROUTINE covmat_degrade


end program degrade_covmat
