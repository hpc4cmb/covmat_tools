! This code reads madam pixel-pixel inverse covariance matrix and a
! given residual noise map. It then computes the chi squared
! test for the given map using the inverse covariance matrix

PROGRAM chisq_test

  USE healpix_types
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, mapfile, mapfile2, form
  CHARACTER(len=filenamelen) :: maskfile
  INTEGER(i4b) :: npixtot, icol, ierr, nhit, rec_len, iArgument
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), ALLOCATABLE :: covmat(:), map1(:), map2(:)
  REAL(dp) :: dof, chisq
  LOGICAL :: there, twomaps=.FALSE., use_mask=.FALSE., pol_only=.FALSE.

  IF (nArguments() < 2)  THEN
     STOP 'Usage: chisq <covmat> <map> [ --map2 <map2>] [--mask <mask>] [--polonly]'
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
     WRITE (*,*) ' Map == ' // TRIM(mapfile)

     maskfile = ''
     iArgument = 3
     DO
        IF (iArgument > nArguments()) EXIT

        CALL getArgument(3, argument)
        IF (INDEX(argument, '-map2') /= 0) THEN
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           mapfile2 = TRIM(argument)
           INQUIRE(file=mapfile2, exist=there)
           IF (.NOT. there) STOP 'mapfile2 does not exist!'
           twomaps = .TRUE.
           WRITE (*,*) ' Got two maps, reporting chi^2 = 2 m1^T N^-1 m2'
           WRITE (*,*) ' Map2 == ' // TRIM(mapfile2)
        ELSE IF (INDEX(argument, '-mask') /= 0) THEN
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           maskfile = TRIM(argument)
           INQUIRE(file=maskfile, exist=there)
           IF (.NOT. there) STOP 'maskfile does not exist!'
           WRITE (*,*) ' Mask == ' // TRIM(maskfile)
           use_mask = .TRUE.
        ELSE IF (INDEX(argument, '-polonly') /= 0) THEN
           WRITE (*,*) ' Will ignore temperature part'
           pol_only = .TRUE.
        ELSE
           WRITE (*,*) 'Unrecognized command line option: ',TRIM(argument)
           STOP
        END IF

        iArgument = iArgument + 1
     END DO
  END IF

  ! Read in the map
  CALL tic
  CALL get_map(mapfile, npixtot, maskfile, pol_only)
  IF (.NOT. twomaps) THEN
     ALLOCATE(map2(0:npixtot-1), stat=ierr)
     IF (ierr /= 0) STOP 'No room for map2'
  END IF
  CALL toc('read map')

  ! Perform (m^T N^-1 m) by storing only a single column of
  ! noise covariance at a time.
  CALL tic
  ALLOCATE(covmat(0:npixtot-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'
  INQUIRE(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  nhit = 0
  DO icol = 0, npixtot - 1
     READ(covmat_unit, rec=icol+1) covmat
     IF (covmat(icol) /= -1) nhit = nhit + 1
     map2(icol) = DOT_PRODUCT(covmat, map1)
  END DO

  IF (twomaps) THEN
     CALL get_map(mapfile2, npixtot, maskfile, pol_only)
  END IF

  chisq = DOT_PRODUCT(map1, map2)
  CALL toc('compute chi-squared')

  IF (twomaps) chisq = 2*chisq

  form = '(a,10es15.5)'
  WRITE (*,form) ' Total chi squared ==  ', chisq

  WRITE (*,*)
  WRITE (*,'(i8,a,i8,a)') nhit, '/', npixtot, ' hit pixels'
  IF (pol_only) THEN
     dof = nhit - nhit/3
  ELSE
     dof = nhit-1
  END IF
  WRITE (*,form) ' chi squared / dof ==  ', chisq/dof
  WRITE (*,*)
  
  WRITE (*,form) ' Expectance == ', REAL(dof,dp)
  WRITE (*,form) '        std == ', (2*dof)**0.5
  WRITE (*,'(/,a,G15.5,a)') &
       ' Deviation from chi squared expected value : ', &
       REAL((chisq-dof),dp)/(2*dof)**0.5, ' sigma'
  WRITE (*,*)


CONTAINS


  SUBROUTINE get_map(file_map, npix, file_mask, pol_only)

    USE pix_tools, ONLY : convert_ring2nest, nside2npix
    USE udgrade_nr, ONLY : udgrade_nest
    USE fitstools, ONLY : input_map, getsize_fits
    
    CHARACTER(len=filenamelen) :: file_map
    INTEGER(i4b) :: npix
    CHARACTER(len=filenamelen) :: file_mask
    LOGICAL(lgt) :: pol_only

    REAL(dp), ALLOCATABLE :: map(:,:), mask(:,:), masktemp(:,:)
    REAL(dp) :: monopole
    INTEGER(i8b) :: dnpix
    INTEGER(i4b) :: npix_mask, nside_mask, nmaps_mask, imap, i, pix, nmaps
    INTEGER(i4b) :: ordering, nside
    
    dnpix = getsize_fits(file_map, nmaps=nmaps, ordering=ordering, nside=nside)
    npix = int(dnpix, i4b)
    ALLOCATE(map(0:npix-1,nmaps), stat=ierr)
    IF (ierr /= 0) STOP 'No room to input map'

    WRITE (*,*) 'Reading ' // TRIM(file_map)
    CALL input_map(file_map, map, npix, nmaps)
    IF (ordering == 1) THEN
       ! Covariance matrix is in NESTED scheme, so must the map be
       WRITE (*,*) ' Note: converting map into NESTED scheme'
       CALL convert_ring2nest(nside, map)
    END IF

    IF (pol_only) THEN
       map(:,1) = 0
    ELSE
       ! Remove T monopole
       monopole = SUM(map(:,1))/npix
       WRITE (*,'(/,a,ES15.5)') '  Removing temperature monopole == ', monopole
       map(:,1) = map(:,1) - monopole
    END IF

    ! FIXME
    ! finish reading the mask and treat properly different nmaps combinations

    ! use a mask?
    IF (LEN_TRIM(file_mask) > 0) THEN
       dnpix = getsize_fits(file_mask, nmaps=nmaps_mask, ordering=ordering, &
            nside=nside_mask)
       npix_mask = int(dnpix, i4b)
       IF (nmaps_mask > nmaps) nmaps_mask = nmaps
       ALLOCATE(masktemp(0:npix_mask-1,nmaps_mask), mask(0:npix-1,nmaps), &
            stat=ierr)
       IF (ierr /= 0) STOP 'No room for mask'
       
       WRITE (*,*) 'Reading ' // TRIM(file_mask)
       CALL input_map(file_mask, masktemp, npix_mask, nmaps_mask)
       IF (ordering == 1) THEN
          ! Covariance matrix is in NESTED scheme, so must the mask be
          WRITE (*,*) ' Note: converting mask into NESTED scheme'
          CALL convert_ring2nest(nside_mask, masktemp)
       END IF

       ! adjust resolution and handle cases where mask has different
       ! number of columns than the map.
       CALL udgrade_nest(masktemp(:,1:nmaps_mask), nside_mask, &
            mask(:,1:nmaps_mask), nside, fmissval=0.0_dp, pessimistic=.TRUE.)
       DO imap = nmaps_mask+1,nmaps
          mask(:,imap) = mask(:,imap-1)
       END DO
       DEALLOCATE(masktemp)

       npix = COUNT(mask == 1)
       WRITE (*,*) ' Mask leaves a total of ',npix,' pixels'
    ELSE
       npix = nmaps * npix
    ENDIF
    
    ! Store the read map into a vector
    ALLOCATE(map1(0:npix-1), stat=ierr)
    IF (ierr /= 0) STOP 'No room for map'
    i = 0
    pix = -1
    imap = 1
    DO
       pix = pix + 1
       IF (pix == nside2npix(nside)) THEN
          pix = 0
          imap = imap + 1
       END IF

       IF (LEN_TRIM(file_mask) > 0) THEN
          IF (mask(pix, imap) == 0) CYCLE
       END IF

       map1(i) = map(pix,imap)

       i = i + 1
       IF (i == npix) EXIT
    END DO

    DEALLOCATE(map)

  END SUBROUTINE get_map

       
END PROGRAM chisq_test
