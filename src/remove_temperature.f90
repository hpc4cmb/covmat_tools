! This code reads in a covariance matrix and writes a new one
! without the II, IQ nor IU blocks
!
! October 9th 2008 Reijo Keskitalo
! 

PROGRAM remove_temperature

  USE healpix_types
  USE extension, ONLY   : nArguments, getArgument
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file
  CHARACTER(len=filenamelen) :: outfilename
  INTEGER(I4B) :: nstokes, nstokes2=0, nside
  INTEGER(I8B) :: npix, icol
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:), covmat2(:)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, rec_len

  IF (nArguments() < 3)  THEN
     STOP 'Usage: merge_covmat <covmat> <nside> <nstokes1> [<covmat_out>]'
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
     WRITE (*,'(a,i4)') 'Nstokes  == ', nstokes

     IF (nArguments() == 4) THEN
        CALL getArgument(4, argument)
        outfilename = TRIM(ADJUSTL(argument))
     ELSE
        outfilename = 'polarization_covmat.dat'
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix = 12*nside**2

  ALLOCATE(covmat(0:npix*nstokes-1), covmat2(0:npix*nstokes2-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for covmat'

  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit,   file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  inquire(iolength=rec_len) covmat(npix+1:)
  OPEN(unit=covmat_unit+2, file=TRIM(ADJUSTL(outfilename)), status='replace', &
       form='unformatted', access='direct', recl=rec_len)

  CALL tic
  DO icol = 1, npix*nstokes
     READ (covmat_unit, rec=icol) covmat
     if (icol > npix) then
        WRITE (covmat_unit+2, rec=icol-npix) covmat(npix+1:)
     end if
  END DO
  CALL toc('trim covmat')

  CLOSE(unit=covmat_unit)
  CLOSE(unit=covmat_unit+2)

  DEALLOCATE(covmat,covmat2)

END PROGRAM remove_temperature
