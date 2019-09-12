! This code reads a pixel-pixel covariance matrix and saves it on disc 
! in expanded (block) form, where correlation of different maps
! is presented in separate blocks
!
!
! Revisions:
! 2012-09-28 : nside and nstokes are now checked internally and cannot
!              be provided

PROGRAM covmat2block_mpi

  USE healpix_types
  USE extension, ONLY : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc, get_matrix_size_dp

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfilename
  INTEGER(i4b) :: nside, nstokes, istokes, isignal, imap
  INTEGER(i8b) :: npixtot, npix, icol, id_read
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: buffer(:)
  LOGICAL :: there, ok
  INTEGER :: ierr, id, ntasks
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: filemode, fileinfo, infilehandle, outfilehandle
  INTEGER(MPI_OFFSET_KIND) :: fileoffset

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (id == 0) THEN
     ok = .TRUE.
     IF (nArguments() < 1 .OR. nArguments() > 2)  THEN
        WRITE(*,'(/,a,/)') &
             'Usage: covmat2block_mpi <covmat> <covmat_out>'
        ok = .FALSE.
     ELSE
        CALL getArgument(1, argument)
        covmat_file = TRIM(argument)
        INQUIRE(file=TRIM(covmat_file), exist=there)
        IF (.NOT. there) THEN
           WRITE (*,*) 'covmat_file does not exist!'
           ok = .FALSE.
        END IF
        WRITE (*,*) ' Reading ' // TRIM(covmat_file)

        outfilename = 'block_'//TRIM(ADJUSTL(covmat_file))
        IF (nArguments() > 1) THEN
           CALL getArgument(2, argument)
           outfilename = TRIM(ADJUSTL(argument))
        END IF
        WRITE (*,*) ' Writing to ' // TRIM(outfilename)
     END IF
  END IF

  CALL mpi_bcast(ok, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
  IF (.NOT. ok) THEN
     CALL mpi_finalize(ierr)
     STOP
  END IF
  CALL mpi_bcast(covmat_file, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, &
       ierr)
  CALL mpi_bcast(outfilename, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, &
       ierr)

  call get_matrix_size_dp(covmat_file, nside, nstokes, npix)
  npixtot = npix
  npix = npix / nstokes
  if (npix <= 0) then
     print *,'ERROR: bad matrix dimensions'
     call mpi_abort(mpi_comm_world, -1, ierr)
  end if
  if (id == 0) then
     WRITE (*,'(a,i10)') ' Npix    == ', npix
     WRITE (*,'(a,i4)') ' Nside   == ', nside
     WRITE (*,'(a,i4)') ' Nstokes == ', nstokes
  end if

  ALLOCATE(buffer(0:npixtot-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for buffer'

  ! Read the matrix full column at a time
  CALL mpi_info_create(fileinfo, ierr)
  fileoffset = 0
  CALL mpi_file_open(mpi_comm_world, TRIM(covmat_file), MPI_MODE_RDONLY, &
       fileinfo, infilehandle, ierr)
  CALL mpi_file_set_view(infilehandle, fileoffset, MPI_REAL8, MPI_REAL8, &
       'native', fileinfo, ierr)

  ! Write the converted matrix single component map at a time
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  CALL mpi_file_open(mpi_comm_world, TRIM(outfilename), filemode, &
       fileinfo, outfilehandle, ierr)
  CALL mpi_file_set_view(outfilehandle, fileoffset, MPI_REAL8, MPI_REAL8, &
       'native', fileinfo, ierr)

  IF (id == 0) CALL tic
  imap = 1
  DO istokes = 1, nstokes
     DO icol = istokes, npixtot, nstokes
        id_read = MODULO(icol, int(ntasks, i8b))
        IF (id /= id_read) THEN
           imap = imap + nstokes
           CYCLE
        END IF
        fileoffset = (icol-1)*npixtot
        CALL mpi_file_read_at(infilehandle, fileoffset, buffer, &
             int(npixtot,i4b), MPI_REAL8, status, ierr)
        DO isignal = 0, nstokes-1
           fileoffset = (imap-1)*npix
           CALL mpi_file_write_at(outfilehandle, fileoffset, &
                buffer(isignal:npixtot-1:nstokes), int(npix,i4b), &
                MPI_REAL8, status, ierr)
           imap = imap + 1
        END DO
     END DO
  END DO

  CALL mpi_file_close(infilehandle, ierr)
  CALL mpi_file_close(outfilehandle, ierr)

  IF (id == 0) CALL toc('convert to block')

  call mpi_finalize(ierr)

END PROGRAM covmat2block_mpi
