! This code reads C binary madam pixel-pixel (inverse) covariance matrix
! writes it as a single column fits file
!
! February 25th 2008 RK

PROGRAM covmat2fits

  USE healpix_types
  USE fitstools, only   : write_bintabh
  USE head_fits, only   : add_card
  USE extension, only   : getArgument, nArguments
  USE covmat_util, only : tic, toc

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfile
  INTEGER(I8B) :: npix, icol
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:, :)
  INTEGER(I4B) :: ierr, nstokes, nside
  LOGICAL :: there
  INTEGER, PARAMETER :: nheaderlines = 128, repeat = 2**20
  CHARACTER(len=80) :: header(1:nheaderlines)

  IF (nArguments()<3)  THEN
     STOP 'Usage: covmat2fits <covmat> <nside> <nstokes> [<outfile>]'
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

     if (nArguments() == 4) then
        CALL getArgument(4, argument)
        outfile = TRIM(argument)
     else
        outfile = TRIM(covmat_file)//'.fits'
     end if
     WRITE (*,*) ' Writing to ' // TRIM(outfile)
  END IF

  npix = 12*nside**2
  ALLOCATE(covmat(0:npix*nstokes-1, 1), stat=ierr)
  IF (ierr/=0) STOP 'no room for covmat'

  OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(covmat_file)), status='old',&
       form='unformatted', access='direct', recl=npix*nstokes*8)
 
  header = ''
  CALL add_card(header, 'creator', 'covmat2fits','covmat software')
  CALL add_card(header, 'nmaps', nstokes, 'number of stokes parameters')
  CALL add_card(header, 'ordering', 'NESTED', 'HEALPIX pixelization scheme')
  CALL add_card(header, 'nside', nside, 'HEALPIX resolution parameter')
  CALL add_card(header, 'units', 'K^2(ant)', 'unit of the elements')
  CALL add_card(header,"HISTORY","Input NCVM = "//TRIM(covmat_file))
  
  call tic
  do icol = 0, nstokes*npix-1
     READ(covmat_unit, rec=icol+1) covmat(:,1)
     CALL write_bintabh(covmat, npix*nstokes, 1, header, nheaderlines, &
          outfile, firstpix=icol*npix*nstokes)!, repeat=repeat)
  end do
  close(covmat_unit)
  call toc('convert matrix')

END PROGRAM covmat2fits
