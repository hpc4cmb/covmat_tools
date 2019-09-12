! This code reads madam pixel-pixel inverse covariance matrix
! and inverts it. The result is stored in both binary and ASCII format

PROGRAM covmat_to_ascii

  USE planck_config, ONLY : iargc, getarg
  USE inputparam
  USE commonparam
  USE output
  USE covmat_util

  IMPLICIT NONE

  EXTERNAL iargc

  CHARACTER(len=80)      :: argument, covmat_file, outfile
  INTEGER(dp)            :: pixel, elements, npix, nstokes, irow, nside
  INTEGER, PARAMETER     :: covmat_unit=55
  INTEGER(dp), PARAMETER :: elements_max=1E8
  REAL(dp), PARAMETER    :: not_read=-1E30
  REAL(dp)               :: row(elements_max)
  REAL(dp), POINTER      :: covmat(:,:)
  INTEGER, POINTER       :: hit_pixels(:)
  LOGICAL                :: there, reduce=.true.
  INTEGER                :: ierr

  IF (iargc/=3)  THEN
     STOP 'Usage: invert_covmat <covmat> <nside> <nstokes>'
  ELSE
     CALL getarg(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP TRIM(covmat_file) // ' does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getarg(2, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getarg(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes
  END IF

  npix = 12*nside**2
  ALLOCATE(covmat(0:npix*nstokes-1, 0:npix*nstokes-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for covmat'

  !OPEN(unit=covmat_unit,file=TRIM(covmat_file),status='old',form='unformatted')

  call tic
  !READ (covmat_unit) covmat
  CALL read_c_matrix(trim(adjustl(covmat_file)), covmat, (nstokes*npix)**2)
  call toc('Read matrix')

  CLOSE(unit=covmat_unit)

  outfile = 'ascii_' // ftrim(covmat_file)

  CALL save_matrix(covmat, filename=trim(adjustl(outfile)), unit=covmat_unit)

  deallocate(covmat)



END PROGRAM covmat_to_ascii
