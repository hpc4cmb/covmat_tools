! This code reads in the eigenvectors and eigenvalues
! and composes the matrix they correspond to
 
PROGRAM evecs2covmat

  USE healpix_types
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, evec_file, eval_file, outfilename
  INTEGER(i4b) :: nstokes, nside
  INTEGER(i8b) :: npix, icol, npix2
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:, :), evec(:, :), evecT(:, :), eval(:)
  REAL(dp) :: limit=1E-6, evalmax
  LOGICAL :: there, invert
  INTEGER(i4b) :: ierr, iArgument

  IF (nArguments() < 4)  THEN
     WRITE (*,'(/,a,/)') 'Usage: evecs2sqrtcovmat <evals> <evecs> ' // &
          '<nside> <nstokes> [--inv] [-o <outfile>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     eval_file = TRIM(argument)
     INQUIRE(file=TRIM(eval_file), exist=there)
     IF (.NOT. there) STOP 'eval_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(eval_file)

     CALL getArgument(2, argument)
     evec_file = TRIM(argument)
     INQUIRE(file=TRIM(evec_file), exist=there)
     IF (.NOT. there) STOP 'evec_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(evec_file)

     CALL getArgument(3, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(4, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     outfilename = 'sqrt_covmat.dat'
     iArgument   = 5
     DO
        CALL getArgument(iArgument, argument)
        IF (INDEX(argument, '-o') /= 0) THEN
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           outfilename = TRIM(ADJUSTL(argument))
        ELSE IF (INDEX(argument, '-inv') /= 0) THEN
           invert = .TRUE.
           WRITE (*,*) 'Inverting while composing'
        ELSE
           WRITE (*,*) 'Unrecognized command line option: ',TRIM(argument)
           STOP
        END IF

        iArgument = iArgument + 1
        IF (iArgument > nArguments()) EXIT
     END DO

     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix = 12*nside**2
  npix2 = npix

  ! read, multiply and write. one column at a time
  CALL tic
  !ALLOCATE(covmat(0:npix*nstokes-1,0:npix*nstokes-1), &
  !     eval(0:npix*nstokes-1), evec(0:npix*nstokes-1,0:npix*nstokes-1), &
  !     evecT(0:npix*nstokes-1,0:npix*nstokes-1), stat=ierr)
  ALLOCATE(covmat(0:npix*nstokes-1,0:npix*nstokes-1), &
       eval(0:npix2*nstokes-1), evec(0:npix*nstokes-1,0:npix2*nstokes-1), &
       evecT(0:npix2*nstokes-1,0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'

  OPEN(unit=covmat_unit, file=TRIM(eval_file), status='old', form='formatted')
  READ(covmat_unit,*) eval
  CLOSE(covmat_unit)
  evalmax = MAXVAL(ABS(eval))

  OPEN(unit=covmat_unit+1, file=TRIM(evec_file), status='old', &
       form='unformatted', access='direct', recl=npix*npix2*nstokes**2*8)
  !     form='unformatted', access='direct', recl=npix**2*nstokes**2*8)
  READ(covmat_unit+1, rec=1) evec
  CLOSE(covmat_unit+1)

  OPEN(unit=covmat_unit+2, file=TRIM(outfilename), status='replace', &
       form='unformatted', access='direct', recl=npix**2*nstokes**2*8)

  covmat    = 0.0
  evec(:,0) = 0.0
  evecT  = TRANSPOSE(evec)
  !write (*,*) 'Omitting first eigenmode'
  !DO icol = 1, npix*nstokes - 1
  DO icol = 0, npix2*nstokes - 1
     IF (ABS(eval(icol))/evalmax < limit) THEN
        evecT(icol,:) = 0.0
        evec(:,icol)  = 0.0
        WRITE (*,'(a,i0)') ' Skipping eigenmode # ',icol
     ELSE
        if (invert) then
           evec(:, icol) = SQRT(1.0/eval(icol))*evec(:, icol)
        else
           evec(:, icol) = SQRT(eval(icol))*evec(:, icol)
        end if
     END IF
  END DO
  covmat = MATMUL(evec, evecT)
  WRITE(covmat_unit+2, rec=1) covmat

  CLOSE(covmat_unit+2)
  CALL toc('compose_covmat')

END PROGRAM evecs2covmat
