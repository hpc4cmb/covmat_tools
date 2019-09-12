! This code reads in the Nobs matrices stored by madam
! and writes an estimate of the covariance matrix based
! purely on white noise characteristics
!
! February 4th 2008 RK

PROGRAM cc2covmat
  use healpix_types
  use extension, only   : getArgument, nArguments
  use fitstools, only   : input_tod
  USE covmat_util, only : tic, toc
  use fitstools
  use head_fits
  use pix_tools

  IMPLICIT NONE

  ! In Madam notation cc means the P^T N^-1_t P matrix
  CHARACTER(len=filenamelen) :: argument, cc_file, covmat_file
  character(len=80) :: header(1:100)
  INTEGER(I8B) :: npix, ipix, icol, isig
  INTEGER(I4B) :: nstokes, nside, ncols, rec_len, ordering
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:, :), covmat_inv(:, :), col(:)
  REAL(dp) :: a, b, c, d, e, f, norm
  LOGICAL :: there, inverse=.false.
  INTEGER(I4B) :: ierr, iarg

  IF (nArguments() < 1)  THEN
     STOP 'Usage: cc2covmat <cc.fits> [-inv] [-o <covmat.dat>]'
  ELSE
     CALL getArgument(1, argument)
     cc_file = TRIM(argument)
     INQUIRE(file=TRIM(cc_file), exist=there)
     IF (.NOT. there) THEN
        WRITE (*,*) TRIM(cc_file) // ' does not exist!'
        STOP
     END IF
     WRITE (*,*) ' Reading ' // TRIM(cc_file)

     iarg = 2
     covmat_file = trim(adjustl(cc_file))//'_covmat.dat'
     do
        if (iarg > nArguments()) exit
        CALL getArgument(iarg, argument)
        if (index(argument,'-o') /= 0) then
           iarg = iarg + 1
           CALL getArgument(iarg, argument)
           covmat_file = TRIM(ADJUSTL(argument))
        else if (index(argument,'-inv') /= 0) then
           inverse = .true.
           write (*,*) ' Computing INVERSE covariance'
        END IF
        iarg = iarg + 1
     end do

     WRITE (*,*) ' Writing to ' // TRIM(covmat_file)
  END IF

  npix = getsize_fits(cc_file, nmaps=ncols, ordering=ordering, nside=nside)

  nstokes = 1
  if (ncols > 1) nstokes = 3

  WRITE (*,*) ' ncols ', ncols
  WRITE (*,*) ' nside ', nside
  WRITE (*,*) ' npix ', npix
  WRITE (*,*) ' ordering ', ordering

  !ncols = 6**((nstokes-1)/2)

  npix = 12*nside**2
  ALLOCATE(covmat(0:npix-1,ncols), covmat_inv(0:npix-1,ncols), &
       col(0:npix*nstokes-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for covmat'

  call tic
  call input_tod(cc_file, covmat_inv, npix, ncols, header)
  call toc('read cc')

  if (ordering == 1) then
     write (*,'(a)') ' Converting map to NESTED'
     call convert_ring2nest(nside, covmat_inv)
  end if

  if (.not. inverse) then
     ! Invert the 3x3 matrices and the whole covmat in the process
     call tic
     do ipix = 0,npix-1
        if (nstokes == 1) then
           covmat(ipix,1) = 1.0/covmat_inv(ipix,1)
        else
           a = covmat_inv(ipix,1)
           b = covmat_inv(ipix,2)
           c = covmat_inv(ipix,3)
           d = covmat_inv(ipix,4)
           e = covmat_inv(ipix,5)
           f = covmat_inv(ipix,6)
           norm = 2.0*b*c*e - c**2*d - b**2*f + a*d*f - a*e**2
           
           covmat(ipix, 1) = (d*f - e**2)/norm
           covmat(ipix, 2) = (c*e - b*f )/norm
           covmat(ipix, 3) = (b*e - c*d )/norm
           covmat(ipix, 4) = (a*f - c**2)/norm
           covmat(ipix, 5) = (c*b - a*e )/norm
           covmat(ipix, 6) = (a*d - b**2)/norm
        end if
     end do
     call toc('invert matrix')
  else
     covmat = covmat_inv
  end if


  ! Finally write the covariance matrix into file
  call tic
  inquire(iolength=rec_len) col
  open (unit=covmat_unit,file=TRIM(covmat_file),status='replace',&
       form='unformatted', access='direct', recl=rec_len)
     do icol = 0, nstokes*npix-1
        col = 0
        ipix = modulo(icol, npix)
        isig = icol / npix + 1
        select case(isig)
        case (1)
           col(ipix)        = covmat(ipix,1)
           col(ipix+  npix) = covmat(ipix,2)
           col(ipix+2*npix) = covmat(ipix,3)
        case (2)
           col(ipix)        = covmat(ipix,2)
           col(ipix+  npix) = covmat(ipix,4)
           col(ipix+2*npix) = covmat(ipix,5)
        case (3)
           col(ipix)        = covmat(ipix,3)
           col(ipix+  npix) = covmat(ipix,5)
           col(ipix+2*npix) = covmat(ipix,6)
        end select
        !col = 0; col(icol)=1+sin(real(ipix)/npix)+1E-3 ! for unit matrix
        write (covmat_unit, rec=icol+1) col
     end do
     close(unit=covmat_unit)
  call toc('Write covmat')

  stop

END PROGRAM cc2covmat
