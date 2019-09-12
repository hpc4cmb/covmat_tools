! This code reads a pixel-pixel covariance matrix and saves it on disc 
! in the expanded (block) form, where correlation of different maps
! is presented in separate blocks

PROGRAM covmat2block
  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=80) :: argument, covmat_file, outfilename
  INTEGER(i4b) :: nstokes, nside, imap, istokes, isignal
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: buffer(:)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, rec_len
  INTEGER(i8b) :: icol, npix

  IF (nArguments() < 3 .OR. nArguments() > 4)  THEN
     STOP 'Usage: covmat2block <covmat> <nside> <nstokes> [<covmat_out>]'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     outfilename = 'block_'//TRIM(ADJUSTL(covmat_file))
     IF (nArguments() > 3) THEN
        CALL getArgument(4, argument)
        outfilename = TRIM(ADJUSTL(argument))
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix = 12*nside**2
  ALLOCATE(buffer(0:npix*nstokes-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for buffer'

  ! Read the matrix full column at a time
  INQUIRE(iolength=rec_len) buffer
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)

  ! Write the converted matrix single component map at a time
  INQUIRE(iolength=rec_len) buffer(0:npix-1)
  OPEN(unit=covmat_unit+1, file=TRIM(ADJUSTL(outfilename)), &
       status='replace', form='unformatted', access='direct', recl=rec_len)

  CALL tic
  imap = 1
  DO istokes = 1, nstokes
     DO icol = istokes, nstokes*npix, nstokes
        READ (covmat_unit, rec=icol) buffer
        DO isignal = 0, nstokes-1
           WRITE (covmat_unit+1, rec=imap) &
                buffer(isignal:npix*nstokes-1:nstokes)
           imap = imap + 1
        END DO
     END DO
  END DO
  CALL toc('convert to block')

  CLOSE(unit=covmat_unit)
  CLOSE(unit=covmat_unit+1)

END PROGRAM covmat2block
