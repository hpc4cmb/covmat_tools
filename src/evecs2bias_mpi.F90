PROGRAM evecs2bias_mpi

!!$ 14 November, 2007 / Torsti Poutanen
!!$ Modified for cmbtools and large files 4 July 2008 / Reijo Keskitalo
  ! Further revisions (RK):
  !   2018-07-06 : Remove T offset before evaluating
  !   MPI parallelization added
  !   mask input addded
  !   2011-04-17 nside and nmaps are now inferred from the eigenvector file

  USE healpix_types
  USE pix_tools
  USE fitstools
  use head_fits
  USE alm_tools
  USE extension, ONLY : getArgument, nArguments
  USE scalapack_tools
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  INTEGER(i4b) :: nside, lmax, l, iargument, imap
  INTEGER(i8b) :: Npix, NNpix
  real(sp), allocatable :: mask(:, :)
  REAL(dp), ALLOCATABLE :: dmapout(:, :), nl(:, :)
  REAL(dp), ALLOCATABLE :: eigenvalue(:), nltest(:, :), w8ring_TQU(:, :)
  REAL(dp), ALLOCATABLE :: eigenvector(:)
  COMPLEX(dpc), ALLOCATABLE :: alm_TGC_test(:, :, :)
  REAL(dp) :: zbounds(1:2), limit=0
  REAL(sp) :: t0, t3, t4
  INTEGER(i4b) :: nmaps, nlheader, ncl, nPlm, ierr, id, ntasks
  CHARACTER(len=80) :: header(1000)
  CHARACTER(len=filenamelen) :: eigenvec_file, eigenval_file, maskfile=''
  CHARACTER(len=filenamelen) :: nlbias_file, argument
  LOGICAL :: there, noInvert=.FALSE.
  REAL(dp), POINTER :: plm(:, :)
  INTEGER(i8b) :: icol, ncol_proc
  INTEGER(i4b) :: status(MPI_STATUS_SIZE), filemode, fileinfo_in, handle_in
  INTEGER(MPI_OFFSET_KIND) :: fileoffset, myfirst, mylast

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)
  IF (id == 0) WRITE (*,'(/,a,i0,a,/)') &
       ' evecs2bias_mpi started with ',ntasks,' tasks.'

  call parse_arguments

  CALL cpu_time(t0)

  Npix  = 12*nside**2
  NNpix = nmaps*Npix
  lmax  = 4*nside

  zbounds = 0

  if (id == 0) then
     PRINT*,''
     PRINT*,'lmin, lmax, nside, Npix, nmaps = ',0,lmax,nside,Npix,nmaps
  end if

  ALLOCATE(w8ring_TQU(1:2*nside, 1:nmaps))
  w8ring_TQU = 1

  ALLOCATE(eigenvalue(0:NNpix-1))
  ALLOCATE(dmapout(0:Npix-1, 1:nmaps))

  if (len_trim(maskfile) > 0) call input_mask(maskfile, mask, npix, id==0)

  ncl = 1
  if (nmaps == 3) ncl = 6
  ALLOCATE(nl(0:lmax, 1:ncl))
  ALLOCATE(alm_TGC_test(1:nmaps, 0:lmax, 0:lmax))
  ALLOCATE(nltest(0:lmax, 1:ncl))

  ALLOCATE(eigenvector(0:NNpix-1))

  nPlm = nside*(lMax+1)*(2*lMax-lMax+2)
  ALLOCATE(plm(0:nPlm-1, nMaps), stat=ierr)
  if (ierr /= 0) stop 'no room to precompute p_lm'
  CALL plm_gen(nSide, lMax, lMax, plm)

  if (id == 0) then
     PRINT*,''
     PRINT*,'Read eigenvalues from  ',TRIM(ADJUSTL(eigenval_file))
     PRINT*,'Read eigenvectors from ',TRIM(ADJUSTL(eigenvec_file))
     OPEN(3,file=eigenval_file)
     READ(3,*) eigenvalue(0:NNpix-1)
     CLOSE(3)
  end if
  call mpi_bcast(eigenvalue, NNpix, MPI_DOUBLE_PRECISION, 0, &
       mpi_comm_world, ierr)

  where (eigenvalue <= limit*maxval(eigenvalue)) eigenvalue = 0

  IF (.NOT. noInvert) THEN
     if (id == 0) write (*,*) 'Inverting the eigenvalues'
     where (eigenvalue /= 0) eigenvalue = 1.0/eigenvalue
  END IF

  if (id == 0) then
     PRINT*,'eigenvalue(0:2) = ',eigenvalue(0:2)
     PRINT*,'eigenvalue(NNpix-1) = ',eigenvalue(NNpix-1)

     PRINT*,''
     PRINT*,'Compute the noise bias'
     CALL cpu_time(t3)
  end if
  nl = 0.0d0


  ! calculate the data distribution
  ncol_proc = NNpix / ntasks
  IF (ncol_proc*ntasks < NNpix) ncol_proc = ncol_proc + 1
  IF (id == 0) WRITE (*,'(" Processing ",i0," columns per task")') ncol_proc

  myfirst = id * ncol_proc
  mylast = (id + 1) * ncol_proc - 1
  if (mylast > NNpix-1) mylast = NNpix-1

  ! open the mpi-io file view INPUT
  CALL mpi_info_create(fileinfo_in, ierr)
  filemode = MPI_MODE_RDONLY
  CALL mpi_file_open(mpi_comm_world, TRIM(eigenvec_file), filemode, &
       fileinfo_in, handle_in, ierr)
  IF (ierr /= 0) STOP 'Failed to open file_in'

  fileoffset = 0
  CALL mpi_file_set_view(handle_in, fileoffset, MPI_REAL8, &
       MPI_REAL8, 'native', fileinfo_in, ierr)
  IF (ierr /= 0) STOP 'Failed to open view to eigenvec_file'

  ! process the eigenvectors
  DO icol = myfirst, mylast
     IF (id == 0 .AND. MODULO(icol-myfirst, ncol_proc/10) == 1) &
          WRITE (*,'(i4,a)') icol*100/(ncol_proc), '% done'

     ! Read in the eigenvector
     fileoffset = icol * NNpix
     CALL mpi_file_read_at(handle_in, fileoffset, eigenvector, NNpix, &
          MPI_REAL8, status, ierr)

     ! Subtract T mean
     eigenvector(:npix) = eigenvector(:npix) - sum(eigenvector(:npix)) / npix

     ! Apply mask
     if (len_trim(maskfile) > 0) then
        do imap = 0, nmaps-1
           eigenvector(imap*npix:(imap+1)*npix-1) = &
                eigenvector(imap*npix:(imap+1)*npix-1) * mask(:, 1)
        end do
     end if

     dmapout = reshape(eigenvector, (/Npix, int(nmaps, i8b)/))
     CALL convert_nest2ring(nside, dmapout)
     if (nmaps == 1) then
        CALL map2alm(nside, lmax, lmax, dmapout(:, 1), alm_TGC_test, zbounds, &
             w8ring_TQU)
     else
        CALL map2alm(nside, lmax, lmax, dmapout, alm_TGC_test, zbounds, &
             w8ring_TQU,plm)
     end if
     CALL alm2cl(lmax, lmax, alm_TGC_test, nltest)
     nl = nl + nltest*eigenvalue(icol)
  END DO

  nltest = nl
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_file_close(handle_in, ierr)

  call mpi_reduce(nltest, nl, (lmax+1)*ncl, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
       mpi_comm_world, ierr)

  if (id == 0) then
     CALL cpu_time(t4)
     PRINT*,'Time to compute the noise_bias [s] = ',t4-t3

     OPEN(3,file=nlbias_file)
     DO l = 0,lmax
        if (ncl == 6) then
           WRITE(3,'(1I5,6E17.9)') &
                l, nl(l,1), nl(l,2), nl(l,3), nl(l,4), nl(l,5), nl(l,6)
        else
           WRITE(3,'(1I5,6E17.9)') l, nl(l,1)
        end if
     END DO
     CLOSE(3)

     ! Write the noise bias into fits file
     header = ''
     call add_card(header, 'creator', 'evecs2bias', &
          'software creating the cls')
     call add_card(header, 'nmaps', nmaps, 'number of stokes parameters')
     call add_card(header, 'lmax', lmax, 'maximum multipole for cls')
     call add_card(header, 'units', 'K^2(ant)', 'unit of the cls')
     call add_card(header,"HISTORY","Input evals = "//TRIM(eigenval_file))
     call add_card(header,"HISTORY","Input evecs = "//TRIM(eigenvec_file))

     call add_card(header, 'TTYPE1', 'TT', 'label for field   1')
     if (ncl == 6) then
        call add_card(header, 'TTYPE2', 'EE', 'label for field   2')
        call add_card(header, 'TTYPE3', 'BB', 'label for field   3')
        call add_card(header, 'TTYPE4', 'TE', 'label for field   4')
        call add_card(header, 'TTYPE5', 'TB', 'label for field   5')
        call add_card(header, 'TTYPE6', 'EB', 'label for field   6')
     end if

     nlheader = SIZE(header)
     if (index(nlbias_file,'.dat') > 0 &
          .or. index(nlbias_file, '.fits') > 0) then
        nlbias_file(len_trim(nlbias_file)-3:len_trim(nlbias_file)+1) = '.fits'
     else
        nlbias_file = TRIM(nlbias_file) // '.fits'
     end if
     call write_asctab(nl, lmax, ncl, header, nlheader, '!'//nlbias_file)
  end if

  call mpi_barrier(mpi_comm_world, ierr)

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)


contains


  subroutine parse_arguments()

    logical :: ok

    integer(i8b) :: npix, nelem

    if (id == 0) then
       ok = .true.
       IF (nArguments() < 3)  THEN
          WRITE (*,'(/,a,/)') &
               'Usage: evecs2bias <eigenval_file> <eigenvec_file> ' // &
               '<nlbias_file> [--noinv] [--limit <limit>] ' // &
               '[--mask <maskfile>]'
          ok = .false.
       ELSE
          CALL getArgument(1, argument)
          eigenval_file = TRIM(argument)
          INQUIRE(file=TRIM(eigenval_file), exist=there)
          IF (.NOT. there) THEN
             WRITE(*,*) TRIM(eigenval_file) // ' does not exist!'
             STOP
          END IF
          WRITE (*,*) ' Reading ' // TRIM(eigenval_file)

          CALL getArgument(2, argument)
          eigenvec_file = TRIM(argument)
          INQUIRE(file=TRIM(eigenvec_file), exist=there)
          IF (.NOT. there) then
             WRITE (*,*) TRIM(eigenvec_file) // ' does not exist!'
             STOP
          END IF
          WRITE (*,*) ' Reading ' // TRIM(eigenvec_file)

          CALL getArgument(3, argument)
          nlbias_file = TRIM(argument)
          WRITE (*,*) ' Writing to ' // TRIM(nlbias_file)

          IF (nArguments() > 3) THEN
             iArgument = 4
             do
                CALL getArgument(iargument, argument)
                IF (INDEX(argument, '-no') /= 0) THEN
                   noInvert = .TRUE.
                   WRITE (*,*) 'Using the eigenvalues without inverse'
                else if (index(argument, '-mask') /= 0) then
                   iArgument = iArgument + 1
                   CALL getArgument(iArgument, argument)
                   maskfile = TRIM(ADJUSTL(argument))
                   inquire(file=trim(maskfile), exist=there)
                   if (.not. there) stop 'maskfile does not exist'
                   write (*,*) 'Applying mask from ',trim(maskfile)
                ELSE IF (INDEX(argument, '-lim') /= 0) THEN
                   iArgument = iArgument + 1
                   CALL getArgument(iargument, argument)
                   READ (argument,*) limit
                   WRITE (*,'("Excluding eigenvalues below ",g20.10)') limit
                ELSE
                   WRITE (*,*) 'Unrecognized argument: '//TRIM(argument)
                   STOP
                END IF
                iArgument = iArgument + 1
                if (iArgument > nArguments()) exit
             end do
             iArgument = iArgument + 1
          END IF
       END IF
    END if

    call mpi_bcast(ok, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    if (.not. ok) then
       call mpi_finalize(ierr)
       stop
    end if
    CALL mpi_bcast(eigenval_file, filenamelen, MPI_CHARACTER, 0, &
         mpi_comm_world, ierr)
    CALL mpi_bcast(eigenvec_file, filenamelen, MPI_CHARACTER, 0, &
         mpi_comm_world, ierr)
    call get_matrix_size_dp(eigenvec_file, nside, nmaps, npix, nelem)
    if (id == 0) then
       WRITE (*,'(a,i4)') ' Nside == ', nside
       WRITE (*,'(a,i4)') ' Nmaps == ', nmaps
       WRITE (*,'(a,i8)') ' Npix == ', npix
       WRITE (*,'(a,i12)') ' Nelem == ', nelem
    end if
    CALL mpi_bcast(nlbias_file, filenamelen, MPI_CHARACTER, 0, &
         mpi_comm_world, ierr)
    CALL mpi_bcast(maskfile, filenamelen, MPI_CHARACTER, 0, &
         mpi_comm_world, ierr)
    call mpi_bcast(noInvert, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    call mpi_bcast(limit, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
  end subroutine parse_arguments


  subroutine input_mask(maskfile, mask, npix, verbose)

    use udgrade_nr, only : udgrade_nest

    character(len=*) :: maskfile
    real(sp), allocatable :: mask(:,:)
    integer(i8b) :: npix
    logical(lgt) :: verbose

    integer(i4b) :: nmaps_mask, ordering_mask, ierr, nside_mask, nside
    integer(i8b) :: npix_mask
    real(sp), allocatable :: masktemp(:,:)
    logical(lgt) :: there

    npix_mask = int(getsize_fits(maskfile, nmaps_mask, ordering_mask), i4b)
    if (verbose) then
       inquire(file=maskfile, exist=there)
       if (.not. there) stop 'maskfile does not exist!'
       write (*,'(a,i0)') ' npix_mask == ', npix_mask
       write (*,'(a,i0)') ' nmaps_mask == ', nmaps_mask
       select case (ordering_mask)
       case (0)
          write (*,'(a)') ' ordering_mask == UNKNOWN'
       case (1)
          write (*,'(a)') ' ordering_mask == RING'
       case (2)
          write (*,'(a)') ' ordering_mask == NEST'
       end select
    end if

    nmaps_mask = 1 ! For now, just use the first component

    allocate(masktemp(0:npix_mask-1, 1:nmaps_mask), stat=ierr)
    if (ierr /= 0) stop 'No room for masktemp'

    allocate(mask(0:npix-1, nmaps_mask), stat=ierr)
    if (ierr /= 0) stop 'No room for mask'

    nside_mask = npix2nside(npix_mask)
    call input_map(maskfile, masktemp, npix_mask, nmaps_mask)
    if (ordering_mask == 1) call convert_ring2nest(nside_mask, masktemp)

    nside = npix2nside(npix)
    call udgrade_nest(masktemp, nside_mask, mask, nside)

    deallocate(masktemp)
  end subroutine input_mask

END PROGRAM evecs2bias_mpi
