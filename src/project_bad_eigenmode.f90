! Compute
! (LNL^t)^{-1} - 
! (LNL^t)^{-1} (Lv) [a^{-2} + (Lv)^t (L N L^t)^{-1} Lv ]^(-1) (Lv)^t (L N L^t)^{-1},
! where
!  N is the noise covariance matrix
!  L is a linear operator, typically some smoothing matrix
!  v is an ill-conditioned eigenvector that we want to project
!    out from the matrix
!
! The resulting matrix can be used to compute chi^2 values that do
! not have contribution from the ill-conditioned eigenmode
!
! 2009-05-18 RK

PROGRAM project_out_bad_eigenmode

  USE healpix_types
  USE fitstools, ONLY   : input_map, getsize_fits
  USE pix_tools, ONLY   : convert_ring2nest
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, mapfile
  INTEGER(i4b) :: nstokes, nside, ordering
  INTEGER(i8b) :: npix, icol
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:), map(:,:), map1(:), map2(:)
  REAL(dp) :: chisq
  LOGICAL :: there
  INTEGER(i4b) :: ierr, isig1

  IF (nArguments() < 2)  THEN
     STOP 'Usage: get_reg_element <inv_covmat> <bad_eigenmode>'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'inv_covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     mapfile = TRIM(argument)
     INQUIRE(file=mapfile, exist=there)
     IF (.NOT. there) STOP 'mapfile does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(mapfile)
  END IF

  ! Read in the map
  CALL tic
  npix = getsize_fits(mapfile, nmaps=nstokes, ordering=ordering, &
       nside=nside)
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
  map2 = map1
  
  CALL tic
  ALLOCATE(covmat(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=npix*nstokes*8)
  DO icol = 0, npix*nstokes - 1
     READ(covmat_unit, rec=icol+1) covmat
     map2(icol) = DOT_PRODUCT(covmat, map1)
  END DO

  chisq = DOT_PRODUCT(map1, map2)
  CALL toc('compute chi-squared')

  write (*,'(/,a,ES15.5,/)') &
       'Without projection, the bad eigenmode has chi^2 = ', chisq

  ! now just compute the regularizing element and subtract it from
  ! the inverse NCM
  map2 = map2 / sqrt(chisq)
  OPEN(unit=covmat_unit+1, file=TRIM(covmat_file)//'.reg', status='replace', &
       form='unformatted', access='direct', recl=npix*nstokes*8)

  call tic
  DO icol = 0, npix*nstokes - 1
     READ(covmat_unit, rec=icol+1) covmat
     map1 = map2*map2(icol)
     WRITE(covmat_unit+1, rec=icol+1) covmat - map1
  END DO
  call toc('project mode out')

  close(unit=covmat_unit)
  close(unit=covmat_unit+1)
       
END PROGRAM project_out_bad_eigenmode
