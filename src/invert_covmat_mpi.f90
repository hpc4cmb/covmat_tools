! Simple interface to SCALAPACK eigenvalue decomposition
! based on the madtools C language interface.
!
! 2008-10-27 Reijo Keskitalo
!
! Revisions:
! 2018-07-08 : Project out the TT offset mode
! 2011-04-17 : nside and nstokes are now checked internally and cannot
!              be provided

PROGRAM invert_covmat_mpi

  USE healpix_types
  use pix_tools, only : nside2npix
  USE extension, ONLY : getArgument, nArguments
  USE covmat_util, only : tic, toc
  use scalapack_tools

  IMPLICIT NONE

  type(scalapack_env) :: spenv
  type(scalapack_matrix_dp) :: evecs, covmat

  CHARACTER(len=filenamelen) :: argument, file_out, file_covmat
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER  :: evals(:)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, id, ntasks, blocksize, nstokes, nside, iargument
  real(dp) :: limit=1E-30
  ! pblas variables
  integer(i4b) :: info
  integer(i8b) :: npix

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  call tic(555)

  IF (nArguments() < 1)  THEN
     IF (id == 0) WRITE(*,*) &
          'Usage: invert_covmat_mpi <file_covmat> '// &
          '[-o <outfile>] [-lim <limit>] [-b <blocksize>]'
     CALL mpi_finalize(ierr)
     STOP
  ELSE IF (id == 0) THEN
     CALL getArgument(1, argument)
     file_covmat = TRIM(argument)
     INQUIRE(file=TRIM(file_covmat), exist=there)
     IF (.NOT. there) THEN
        write (*,'(a)') 'file_covmat does not exist: ' // trim(file_covmat)
        call mpi_abort(mpi_comm_world, -1, ierr)
     END IF
     WRITE (*,'(/,a)') ' Reading matrix from ' // TRIM(file_covmat)

     file_out = trim(file_covmat) // '.inverse'
     blocksize = 32
     iArgument = 2
     do
        if (iArgument > nArguments()) exit

        CALL getArgument(iArgument, argument)
        if (index(argument, '-o') /= 0) then
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           file_out = TRIM(ADJUSTL(argument))
        else if (index(argument, '-lim') /= 0) then
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           read(argument, *, iostat=ierr) limit
           if (ierr /= 0) then
              print *,'Failed to parse limit from ' // trim(argument)
              stop
           end if
           write (*,*) 'Using threshold eval > evalmax * ',limit
        else if (index(argument, '-b') /= 0) then
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           read(argument, *, iostat=ierr) blocksize
           if (ierr /= 0) then
              print *,'Failed to parse blocksize from ' // trim(argument)
              stop
           end if
        else
           write (*,*) 'Unrecognized command line option: ', trim(argument)
           stop
        end if

        iArgument = iArgument + 1
     end do

     write (*,'(" Block == ",i4)') blocksize
     WRITE (*,'(a,/)') ' Writing to ' // TRIM(file_out)

  END IF

  CALL mpi_bcast(file_covmat, filenamelen, mpi_character, 0, mpi_comm_world, &
       ierr)
  call get_matrix_size_dp(file_covmat, nside, nstokes, npix)
  if (npix <= 0) then
     print *,'ERROR: bad matrix dimensions'
     call mpi_abort(mpi_comm_world, -1, ierr)
  end if
  if (id == 0) then
     WRITE (*,'(a,i10)') ' Npix    == ', npix
     WRITE (*,'(a,i4)') ' Nside   == ', nside
     WRITE (*,'(a,i4)') ' Nstokes == ', nstokes
  end if
  CALL mpi_bcast(blocksize, 1, mpi_integer, 0 ,mpi_comm_world, ierr)
  CALL mpi_bcast(limit, 1, mpi_double_precision, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(file_out, filenamelen,mpi_character, 0, mpi_comm_world, ierr)

  call init_scalapack(spenv, blocksize)

  info = 0

  ! read in the matrix
  if (ID == 0) call tic
  call read_scalapack_matrix_dp(spenv, covmat, file_covmat, npix, npix)
  CALL mpi_barrier(mpi_comm_world, ierr)
  if (ID == 0) call toc('Read matrix')

  ! symmetrize the matrix
  if (ID == 0) call tic
  call symmetrize_scalapack_matrix_dp(spenv, covmat)
  CALL mpi_barrier(mpi_comm_world, ierr)
  if (ID == 0) call toc('Symmetrize matrix')

  ! decompose the matrix
  if (ID == 0) call tic
  call eigendecomp_scalapack_matrix_dp(spenv, covmat, evals, evecs)
  CALL mpi_barrier(mpi_comm_world, ierr)
  if (ID == 0) call toc('Eigen-decompose matrix')

  call project_tt_offset()

  ! save the decomposition
  if (id == 0) then
     open(file=trim(file_covmat)//'.eigenvals', form='formatted', &
          status='replace', unit=55)
     write(55,'(es20.10)') evals
     close(55)
  end if
  if (id == 0) call tic
  call write_scalapack_matrix_dp(spenv, evecs, trim(file_covmat)//'-evecs')
  if (ID == 0) call toc('Save eigenvectors')


  ! invert the eigenvalues for recomposing
  where(evals < maxval(evals)*limit) evals = 0
  if (id == 0) then
     open(file=trim(file_covmat)//'.eigen', form='formatted', &
          status='replace', unit=55)
     write(55,'(es20.10)') evals
     close(55)
  end if
  where(evals /= 0) evals = 1.0 / sqrt(evals)

  call apply_eigenvalues()

  ! Multiply
  if (ID == 0) call tic
  call pdgemm('No transpose', 'Transpose', npix, npix, npix, &
       1.0_dp, evecs%data, 1, 1, evecs%desc, &
       evecs%data, 1, 1, evecs%desc, &
       0.0_dp, covmat%data, 1, 1, covmat%desc)
  if (ID == 0) call toc('Multiplication')


  if (ID == 0) call tic
  call write_scalapack_matrix_dp(spenv, covmat, file_out)
  if (ID == 0) call toc('Save inverse')

  if (id == 0) call toc('invert_covmat_mpi',555)


  CALL mpi_finalize(ierr)


contains


  subroutine apply_eigenvalues()
    ! apply the eigenvalues by looping over the local array
    ! and converting locations into global column numbers
    integer(i4b) :: ielement, icolumn_tot, irow, icol
    if (ID == 0) call tic
    ielement = 0
    do icol = 0, evecs%ncol - 1
       icolumn_tot = &
            (int(icol/spenv%block)*spenv%npcol + spenv%mycol)*spenv%block &
            + modulo(icol, spenv%block)
       do irow = 0, evecs%nrow - 1
          ielement = ielement + 1
          evecs%data(ielement) = evecs%data(ielement) * evals(icolumn_tot+1)
       end do
    end do
    if (ID == 0) call toc('Apply eigenvalues')
  end subroutine apply_eigenvalues


  subroutine project_tt_offset()
    ! Project the TT offset out from each eigenvector
    integer(i4b) :: ielement, icolumn_tot, irow_tot, irow, icol
    real(dp), allocatable :: evec_mean(:)
    if (ID == 0) call tic
    allocate(evec_mean(0:npix*nstokes-1), stat=ierr)
    if (ierr /= 0) then
       print *, 'No room for evec_mean'
       call mpi_abort(mpi_comm_world, -1, ierr)       
    end if
    ! First sum the T elements in each eigenvector
    evec_mean = 0
    ielement = 0
    do icol = 0, evecs%ncol - 1
       icolumn_tot = &
            (int(icol/spenv%block)*spenv%npcol + spenv%mycol)*spenv%block &
            + modulo(icol, spenv%block)
       do irow = 0, evecs%nrow - 1
          irow_tot = &
               (int(irow/spenv%block)*spenv%nprow + spenv%myrow)*spenv%block &
               + modulo(irow, spenv%block)
          ielement = ielement + 1
          ! Skip diagonal elements in measuring the T mean
          if (irow_tot < npix .and. irow_tot /= icolumn_tot) then
             evec_mean(icolumn_tot) = evec_mean(icolumn_tot) &
                  + evecs%data(ielement)
          end if
       end do
    end do
    call mpi_allreduce(MPI_IN_PLACE, evec_mean, npix*nstokes, &
         MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
    ! ... then average ...
    evec_mean = evec_mean / (npix - 1)
    ! ... and finally subtract the T mean from each vector
    ielement = 0
    do icol = 0, evecs%ncol - 1
       icolumn_tot = &
            (int(icol/spenv%block)*spenv%npcol + spenv%mycol)*spenv%block &
            + modulo(icol, spenv%block)
       do irow = 0, evecs%nrow - 1
          irow_tot = &
               (int(irow/spenv%block)*spenv%nprow + spenv%myrow)*spenv%block &
               + modulo(irow, spenv%block)
          ielement = ielement + 1
          if (irow_tot < npix) then
             evecs%data(ielement) = evecs%data(ielement) &
                  - evec_mean(icolumn_tot)
          end if
       end do
    end do
    deallocate(evec_mean)
    if (ID == 0) call toc('Project TT offset')
  end subroutine project_tt_offset



END PROGRAM invert_covmat_mpi
