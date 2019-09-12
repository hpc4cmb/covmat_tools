!> @file
!! Converts a C binary matrix into F90 binary as required
!! by the Commander code

!> Usage: covmat2commander <covmat> <npix_tot> <nstokes> <multiplier> [<outfile>]
!!
!! @param covmat input C flat binary file containing the matrix
!! @param npix_tot total number of I+Q+U pixels
!! @param nstokes number of Stokes parameters for the matrix (1 or 3)
!! @param outfile optional name of the output file (default is <covmat>.commander)
PROGRAM covmat2commander

  USE healpix_types
  USE pix_tools, only : nside2npix
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfilename
  INTEGER(i4b) :: nstokes, rec_len
  INTEGER(i8b) :: icol, npix_tot
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:)
  REAL(dp) :: multiplier
  LOGICAL :: there
  INTEGER(i4b) :: ierr

  IF (nArguments() < 4)  THEN
     WRITE (*,'(/,a,/)') '  Usage: covmat2commander <covmat> <npix_tot> <nstokes> ' // &
          '<multiplier> [<outfile>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *, iostat=ierr) npix_tot
     if (ierr /= 0) stop 'Unable to parse total pixels'
     WRITE (*,'(a,i4)') 'Npix_tot == ', npix_tot

     CALL getArgument(3, argument)
     READ (argument, *, iostat=ierr) nstokes
     if (ierr /= 0) stop 'Unable to parse nstokes'
     WRITE (*,'(a,i4)') 'Nstokes  == ', nstokes

     CALL getArgument(4, argument)
     READ (argument, *, iostat=ierr) multiplier
     if (ierr /= 0) stop 'Unable to parse multiplier'
     WRITE (*,*) ' Multiplying the covariance matrix by ', multiplier

     outfilename = TRIM(ADJUSTL(covmat_file))//'.commander'
     IF (nArguments() > 4) THEN
        CALL getArgument(5, argument)
        outfilename = TRIM(ADJUSTL(argument))
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

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
       form='unformatted')

  ! First the header lines
  WRITE(covmat_unit+1) npix_tot
  ! two other integers that do not seem be used for anything right now..
  WRITE(covmat_unit+1) 2 ! ordering is always NESTED
  IF (nstokes == 1) THEN
     WRITE(covmat_unit+1) 0 ! unpolarized
  ELSE
     WRITE(covmat_unit+1) 1 ! polarized
  END IF

  DO icol = 1, npix_tot
     READ(covmat_unit, rec=icol) covmat
     IF (multiplier /= 1) covmat = covmat * multiplier
     WRITE(covmat_unit+1) covmat
  END DO

  CLOSE(covmat_unit)
  CLOSE(covmat_unit+1)
  CALL toc('covmat2commander')

END PROGRAM covmat2commander
