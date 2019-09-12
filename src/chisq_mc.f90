! This code reads madam pixel-pixel inverse covariance matrix and a
! given residual noise map. It then computes the chi squared
! test for the given map using the inverse covariance matrix

PROGRAM chisq_test

  USE healpix_types
  USE fitstools, ONLY   : input_map, getsize_fits
  USE pix_tools, ONLY   : convert_ring2nest
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments
  USE rngmod

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, mapfile, mapfile2, form
  INTEGER(i4b) :: pixel, npix, nstokes, icol, nside, ordering
  INTEGER(i8b) :: dnpix
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:), map(:, :), map1(:), map2(:)
  REAL(dp) :: chisq, sigma=-1E30
  REAL(dp) :: dof, gauss
  LOGICAL :: there, twomaps=.false.
  INTEGER(i4b) :: ierr, isig1, nhit, seed, rec_len
  TYPE(planck_rng) :: rnghandle

  IF (nArguments() < 2)  THEN
     STOP 'Usage: chisq <covmat> <map> [<sigma>]'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     mapfile = TRIM(argument)
     INQUIRE(file=mapfile, exist=there)
     IF (.NOT. there) STOP 'mapfile does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(mapfile)

     if (nArguments() > 2) then
        CALL getArgument(3, argument)
        read(argument, *) sigma
        WRITE (*,'(a,es15.5)') &
             '  Will add white noise, sigma = ',sigma
     end if
  END IF

  ! Read in the map
  CALL tic
  dnpix = getsize_fits(mapfile, nmaps=nstokes, ordering=ordering, &
       nside=nside)
  npix = int(dnpix, i4b)
  ALLOCATE(map(0:npix-1,nstokes), map1(0:npix*nstokes-1), &
       map2(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for maps'
  CALL input_map(mapfile, map, npix, nstokes)
  IF (ordering == 1) THEN
     ! Covariance matrix is in NESTED scheme, so must the map be
     WRITE (*,*) ' Note: converting map into NESTED scheme'
     CALL convert_ring2nest(nside, map)
  END IF
  CALL toc('read map')

  ! Store the read map into a vector
  do isig1 = 1, nstokes
     map1((isig1-1)*npix:isig1*npix-1) = map(:,isig1)
  end do

  ! make noise.
  if (sigma > 0) then
     CALL SYSTEM_CLOCK(count=seed)
     CALL rand_init(rnghandle, seed)
     DO pixel = 0, npix*nstokes-1
        gauss       = rand_gauss(rnghandle)
        map1(pixel) = map1(pixel) + gauss*sigma
     END DO
  end if
  
  map2 = map1

  ! Perform (m^T N^-1 m) by storing only a single column of
  ! noise covariance at a time. To report a breakdown of the
  ! chisq, this is done in blocks
  CALL tic
  ALLOCATE(covmat(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  nhit = 0
  DO icol = 0, npix*nstokes - 1
     READ(covmat_unit, rec=icol+1) covmat
     IF (icol < npix .AND. covmat(icol) > 0) nhit = nhit + 1
     map2(icol) = DOT_PRODUCT(covmat, map1)
  END DO

  if (twomaps) then
     CALL input_map(mapfile2, map, npix, nstokes)
     IF (ordering == 1) THEN
        ! Covariance matrix is in NESTED scheme, so must the map be
        WRITE (*,*) ' Note: converting map into NESTED scheme'
        CALL convert_ring2nest(nside, map)
     END IF
     do isig1 = 1, nstokes
        map1((isig1-1)*npix:isig1*npix-1) = map(:,isig1)
     end do
  end if

  chisq = DOT_PRODUCT(map1, map2)
  CALL toc('compute chi-squared')

  if (twomaps) chisq = 2*chisq

  form = '(a,10es15.5)'
  WRITE (*,form) ' Total chi squared ==  ', chisq

  WRITE (*,*)
  WRITE (*,'(i8,a)') nhit, ' hit pixels'
  dof = nhit*nstokes-1
  WRITE (*,form) ' chi squared / dof ==  ', chisq/dof
  WRITE (*,*)
  
  WRITE (*,form) ' Expectance == ', REAL(dof,dp)
  WRITE (*,form) '        std == ', (2*dof)**0.5
  WRITE (*,'(/,a,G15.5,a)') &
       ' Deviation from chi squared expected value : ', &
       REAL((chisq-dof),dp)/(2*dof)**0.5, ' sigma'
  WRITE (*,*)
       
END PROGRAM chisq_test
