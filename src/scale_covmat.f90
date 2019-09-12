! This simple code multiplies the covariance matrix
! by a given value
 
PROGRAM scale_covmat

  USE healpix_types
  USE pix_tools, only : nside2npix
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfilename
  INTEGER(i4b) :: nstokes, nside, rec_len
  INTEGER(i8b) :: icol, npix_tot
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:)
  REAL(dp) :: multiplier
  LOGICAL :: there
  INTEGER(i4b) :: ierr

  IF (nArguments() < 4)  THEN
     WRITE (*,'(/,a,/)') '  Usage: scale_covmat <covmat> <nside> <nstokes> ' // &
          '<multiplier> [<outfile>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *, iostat=ierr) nside
     if (ierr /= 0) stop 'Unable to parse nside'
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(3, argument)
     READ (argument, *, iostat=ierr) nstokes
     if (ierr /= 0) stop 'Unable to parse nstokes'
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     CALL getArgument(4, argument)
     READ (argument, *, iostat=ierr) multiplier
     if (ierr /= 0) stop 'Unable to parse multiplier'
     WRITE (*,*) ' Multiplying the covariance matrix by ', multiplier

     outfilename = 'scaled_'//TRIM(ADJUSTL(covmat_file))
     IF (nArguments() > 4) THEN
        CALL getArgument(5, argument)
        outfilename = TRIM(ADJUSTL(argument))
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix_tot = nside2npix(nside) * nstokes

  ! read, multiply and write. one column at a time
  CALL tic
  ALLOCATE(covmat(npix_tot), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'

  ! input
  INQUIRE(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)

  ! output
  OPEN(unit=covmat_unit+1, file=TRIM(outfilename), status='replace', &
       form='unformatted', access='direct', recl=rec_len)

  ! do the scaling
  DO icol = 1, npix_tot
     READ(covmat_unit, rec=icol) covmat
     WRITE(covmat_unit+1, rec=icol) covmat*multiplier
  END DO

  CLOSE(covmat_unit)
  CLOSE(covmat_unit+1)
  CALL toc('scale_covmat')

END PROGRAM scale_covmat
