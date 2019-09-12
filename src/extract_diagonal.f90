!> @ file
!! Dump the covariance matrix diagonal into fits file

!> Usage: Usage: extract_diagonal <covmat> <nside> <nstokes> [<outfile.fits>]
!!
!! @param covmat Input covariance matrix
!! @param nside Covariance matrix resolution
!! @param nstokes Number of Stokes parameters (1 or 3)
!! @param outfile.fits (optional) name of the output fits file. Default: <covmat>_diagonal.fits
PROGRAM extract_diagonal

  USE healpix_types
  USE pix_tools, only : nside2npix
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY : getArgument, nArguments
  USE head_fits, only : write_minimal_header
  use fitstools, only : output_map

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfilename
  INTEGER(i8b) :: npix, npix_tot, ipix, icol
  INTEGER(i4b) :: istokes, nstokes, nside, rec_len
  INTEGER(i4b), PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:), diagonal(:, :)
  LOGICAL(lgt) :: there
  INTEGER(i4b) :: ierr
  character(len=80) :: header(1000)

  IF (nArguments() < 3 .or. nArguments() > 4)  THEN
     WRITE (*,'(/,a,/)') ' Usage: extract_diagonal <covmat> <nside> <nstokes> ' // &
          '[<outfile.fits>]'
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

     outfilename = '!'//TRIM(ADJUSTL(covmat_file))//'_diagonal.fits'
     IF (nArguments() == 4) THEN
        CALL getArgument(4, argument)
        outfilename = TRIM(ADJUSTL(argument))
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix = nside2npix(nside)
  npix_tot = npix * nstokes

  ! read the matrix one column at a time and save the diagonal
  CALL tic
  ALLOCATE(covmat(npix_tot), diagonal(0:npix-1, nstokes), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'

  INQUIRE(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)

  istokes = 1
  icol = 1
  do istokes = 1, nstokes
     do ipix = 0, npix-1
        READ(covmat_unit, rec=icol) covmat
        diagonal(ipix, istokes) = covmat(icol)
        icol = icol + 1
     end do
  END DO

  CLOSE(covmat_unit)

  ! output the results
  header = ' '
  call write_minimal_header(header, 'MAP', nside=nside, ordering='NESTED', &
       creator='extract_diagonal', polar=(nstokes==3))
  call output_map(diagonal, header, outfilename)

  CALL toc('extract diagonal')

END PROGRAM extract_diagonal
