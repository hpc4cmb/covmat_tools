! This code reads madam alm-alm covariance matrix and inverts it. 
! The result is stored in both binary and ASCII format
!
! Matrix is assumed not to include ell < 2, m < 0 modes

PROGRAM invert_harmonic_covmat

  USE healpix_types
  USE extension, ONLY   : nArguments, getArgument
  USE covmat_util, ONLY : tic, toc, eigen_invert_matrix

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfile
  INTEGER(I4B) :: nalm, nalm2, nstokes, lmin, lmax
  INTEGER, PARAMETER :: covmat_unit=55
  COMPLEX(dpc), POINTER :: covmat(:, :)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, icol, rec_len

  IF (nArguments() /= 4)  THEN
     STOP 'Usage: invert_covmat <covmat> <lmin> <lmax> <nstokes>'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) lmin
     WRITE (*,'(a,i4)') '  lMin    == ', lmin

     CALL getArgument(3, argument)
     READ (argument, *) lmax
     WRITE (*,'(a,i4)') '  lMax    == ', lmax

     CALL getArgument(4, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') '  Nstokes == ', nstokes
  END IF

  nalm = (lmax+1)**2 - lmin**2
  nalm2 = nalm - (nalm - (lmax+1-lmin))/2 ! exclude degenerate modes
  ALLOCATE(covmat(nalm2*nstokes, nalm2*nstokes), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat'

  CALL tic
  inquire(iolength=rec_len) covmat(:, 1)
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old',&
       form='unformatted', access='direct', recl=rec_len)
  DO icol = 1, nalm2*nstokes
     READ (covmat_unit, rec=icol) covmat(:, icol)
  END DO
  CLOSE(unit=covmat_unit)
  CALL toc('Read matrix')
  
  CALL tic
  CALL eigen_invert_matrix(covmat)
  CALL toc('Invert matrix')

  CALL tic
  outfile = TRIM(ADJUSTL(TRIM(covmat_file)//'.inverse'))
  WRITE (*,*) ' Writing to: ', TRIM(outfile)
  OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(outfile)), status='replace', &
       form='unformatted', access='direct', recl=rec_len)
  DO icol = 1, nalm2*nstokes
     WRITE (covmat_unit, rec=icol) covmat(:, icol)
  END DO
  CLOSE(unit=covmat_unit)
  CALL toc('Write matrix')

  DEALLOCATE(covmat)

END PROGRAM invert_harmonic_covmat
