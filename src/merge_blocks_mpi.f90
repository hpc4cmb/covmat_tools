! Command line utility to merge the blocks of the covariance matrix
! after the massive parallel run.
!
! input  : fortran binary blocks from madam_NCVM run
! output : a single monolithic inverse covariance matrix file, C binary


PROGRAM merge_blocks_mpi

  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(len=1024) :: argument, covmat_file, first_file, outfilename
  INTEGER(i4b) :: nstokes, nside
  INTEGER(i8b) :: npix
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat_in(:,:), covmat_out(:,:)
  LOGICAL :: there, ok
  INTEGER(i4b) :: ierr, icol_in, icol_out, n, id, ntasks, id_read
  INTEGER(i4b) :: status(MPI_STATUS_SIZE), rec_len

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (id == 0) THEN
     ok = .TRUE.
     IF (nArguments() /= 3)  THEN
        WRITE (*,'(/,a,/)') &
             ' Usage: merge_blocks <first_block> <nside> <nstokes>'
        ok = .FALSE.
     ELSE
        WRITE (*,'(/,a,i0,a,/)') &
             ' merge_blocks started with ',ntasks,' tasks.'
        IF (ntasks < 2) THEN
           WRITE (*,*) ' ERROR, need as least 2 tasks'
           ok = .FALSE.
        END IF
        
        CALL getArgument(1, argument)
        first_file = TRIM(argument)
        INQUIRE(file=TRIM(first_file), exist=there)
        IF (.NOT. there) THEN
           WRITE (*,*) 'covmat_file does not exist!'
           ok = .FALSE.
        END IF
        WRITE (*,*) ' Reading ' // TRIM(first_file)

        CALL getArgument(2, argument)
        READ (argument, *) nside
        WRITE (*,'(a,i4)') 'Nside   == ', nside
        IF (MODULO(nside,2) /= 0) THEN
           WRITE (*,*) 'Illegal nside'
           ok = .FALSE.
        END IF

        CALL getArgument(3, argument)
        READ (argument, *) nstokes
        WRITE (*,'(a,i4)') 'Nstokes == ', nstokes
        IF (.NOT. (nstokes == 1 .OR. nstokes == 3)) THEN
           WRITE (*,*) 'Illegal nstokes'
           ok = .FALSE.
        END IF

        CALL get_outfilename(first_file, outfilename)
        WRITE (*,*) 'Writing to ' // TRIM(outfilename)
     END IF
  END IF

  CALL mpi_bcast(ok, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
  IF (.NOT. ok) THEN
     CALL mpi_finalize(ierr)
     STOP
  END IF
  CALL mpi_bcast(first_file,  1024, MPI_CHARACTER, 0, mpi_comm_world, ierr)
  !CALL mpi_bcast(outfilename, 1024, MPI_CHARACTER, 0, mpi_comm_world, ierr)
  CALL get_outfilename(first_file, outfilename)
  CALL mpi_bcast(nside,   1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(nstokes, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)

  npix = 12*nside**2
  n = nstokes**2

  ALLOCATE(covmat_in(n, 0:npix-1), covmat_out(0:npix*nstokes-1, nstokes), &
       stat=ierr) ! need only one column
  IF (ierr /= 0) STOP 'no room for covmat'

  IF (id == 0) THEN
     CALL tic
     inquire(iolength=rec_len) covmat_out
     OPEN(unit=covmat_unit+1, file=TRIM(ADJUSTL(outfilename)), &
          status='replace', form='unformatted', access='direct', &
          recl=npix*nstokes**2*8)
  END IF

  CALL mpi_barrier(mpi_comm_world, ierr)
  icol_in = 1
  icol_out = 1
  covmat_file = first_file
  loop_files : DO 
     ! cycle through available files reading the same column from
     ! all of the files and writing it into a single file
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (there) THEN
        id_read = MODULO(icol_out-1, ntasks-1)+1
        IF (id == id_read) THEN
           inquire(iolength=rec_len) covmat_in
           OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old',&
                form='unformatted', access='direct', recl=rec_len)
           READ (covmat_unit, rec=icol_in) covmat_in
           CLOSE(unit=covmat_unit)

           IF (nstokes == 1) THEN
              covmat_out(:,1) = covmat_in(1,:)
           ELSE
              covmat_out(0:npix*nstokes-1:nstokes, 1) = covmat_in(1, :)
              covmat_out(1:npix*nstokes-1:nstokes, 1) = covmat_in(2, :)
              covmat_out(2:npix*nstokes-1:nstokes, 1) = covmat_in(3, :)
              
              covmat_out(0:npix*nstokes-1:nstokes, 2) = covmat_in(4, :)
              covmat_out(1:npix*nstokes-1:nstokes, 2) = covmat_in(5, :)
              covmat_out(2:npix*nstokes-1:nstokes, 2) = covmat_in(6, :)

              covmat_out(0:npix*nstokes-1:nstokes, 3) = covmat_in(7, :)
              covmat_out(1:npix*nstokes-1:nstokes, 3) = covmat_in(8, :)
              covmat_out(2:npix*nstokes-1:nstokes, 3) = covmat_in(9, :)
           END IF

           CALL mpi_send(covmat_out, npix*nstokes**2, MPI_DOUBLE_PRECISION, &
                0, icol_out, mpi_comm_world, ierr)
        END IF

        IF (id == 0) THEN
           CALL mpi_recv(covmat_out, npix*nstokes**2, MPI_DOUBLE_PRECISION, &
                id_read, icol_out, mpi_comm_world, status, ierr)
           WRITE (covmat_unit+1, rec=icol_out) covmat_out
        END IF

        icol_out = icol_out + 1
        IF (icol_out > npix) EXIT loop_files

        CALL update_covmat_file(covmat_file)
     ELSE
        covmat_file = first_file
        icol_in     = icol_in + 1
     END IF
  END DO loop_files

  CALL mpi_barrier(mpi_comm_world, ierr)

  IF (id == 0) THEN
     CLOSE(unit=covmat_unit+1)
     CALL toc('merge blocks')
  END IF

  CALL mpi_finalize(ierr)


CONTAINS


  SUBROUTINE get_outfilename(covmat_file, outfilename)
    CHARACTER(len=*) :: covmat_file, outfilename
    INTEGER :: istart, iend

    istart = INDEX(covmat_file, '_block_')
    IF (istart < 0) STOP 'bad block number'
    iend = istart + 7 + INDEX(covmat_file(istart+7:), '_')

    outfilename = ' '
    outfilename = covmat_file(1:istart) // TRIM(covmat_file(iend:))
        
  END SUBROUTINE get_outfilename


  SUBROUTINE update_covmat_file(covmat_file)
    CHARACTER(len=80) :: covmat_file, stemp
    INTEGER           :: istart, iend, blok, ilength

    istart = INDEX(covmat_file, '_block_') + 6  ! index before block number
    IF (istart-7 < 0) STOP 'bad block number'
    iend = INDEX(covmat_file(istart+1:80), '_') ! index afer block number
    ilength = iend-1
    iend = iend + istart

    stemp = covmat_file(istart+1:iend-1)
    READ(stemp, *) blok
    blok = blok + 1
    WRITE (stemp, *) blok

    covmat_file = covmat_file(1:istart) &
         // TRIM(ADJUSTL(stemp)) &
         // covmat_file(iend:80)

  END SUBROUTINE update_covmat_file

END PROGRAM merge_blocks_mpi
