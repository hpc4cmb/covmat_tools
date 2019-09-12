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

  CHARACTER(len=80)  :: argument, covmat_file, outfile
  INTEGER(dp)        :: pixel, npix
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER  :: covmat(:)
  INTEGER            :: ierr, row, col, nstokes, nside
  LOGICAL            :: there

  IF (iargc<3)  THEN
     STOP 'Usage: covmat_to_ascii <covmat> <nside> <nstokes> [<outfile>]'
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

     if (iargc == 4) then
        CALL getarg(4, argument)
        outfile = TRIM(argument)
     else
        outfile = 'ascii_'//covmat_file
     end if
     WRITE (*,*) ' Writing to ' // TRIM(outfile)
  END IF

  npix = 12*nside**2
  ALLOCATE(covmat(0:npix*nstokes-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for covmat'

  OPEN(unit=covmat_unit,  file=TRIM(covmat_file),status='old',form='unformatted')
  open(unit=covmat_unit+1,file=trim(outfile), status='replace')

  call tic
  do col = 0,nstokes*npix-1
     READ (covmat_unit) covmat
     do row = 0, nstokes*npix-1
        WRITE (covmat_unit+1, '(ES15.5)', advance='no') &
             covmat(row)
     end do
     write (covmat_unit+1,*)
  end do

  CLOSE(unit=covmat_unit)
  CLOSE(unit=covmat_unit+1)
  call toc('convert matrix')

  deallocate(covmat)

END PROGRAM covmat_to_ascii
