! This is a simple command line utility to check map statistics
! like the mean and standard deviation. When supplied with
! two maps, it analyzes the difference of the two maps.
!
! 02/02/2010 RK

PROGRAM mapstats

  USE healpix_types
  USE fitstools, ONLY   : input_map, getsize_fits
  USE pix_tools, ONLY   : convert_ring2nest
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, mapfile, mapfile2, maskfile
  INTEGER(i4b) :: nstokes, nside, nside2, nstokes2
  INTEGER(i8b) :: npix, npix2, ipix
  REAL(dp), pointer :: map(:, :)=>null(), map1(:, :)=>null(), map2(:, :)=>null()
  REAL(dp) :: mean, stdev, mapmin, mapmax
  LOGICAL :: there, twomaps=.false.
  LOGICAL, POINTER :: mask(:)=>null()
  INTEGER(i4b) :: ierr, nhit, istokes, iarg
  character(len=1) :: stokes(3) = (/"I", "Q", "U"/)

  IF (nArguments() < 1)  THEN
     STOP 'Usage: mapstats <map> [<map2>] [--mask <mask>]'
  ELSE
     CALL getArgument(1, argument)
     mapfile = TRIM(argument)
     INQUIRE(file=mapfile, exist=there)
     IF (.NOT. there) STOP 'mapfile does not exist!'
     WRITE (*,*) ' Map1 == ' // TRIM(mapfile)

     iarg  = 2
     maskfile = ''
     do
        if (iarg > narguments()) exit
        call getArgument(iarg, argument)
        if (index(argument, '-mask') > 0) then
           iarg = iarg + 1
           call getArgument(iarg, argument)
           maskfile = trim(adjustl(argument))
           inquire(file=maskfile, exist=there)
           if (.not. there) then
              write (*,'(a)') 'ERROR: ' // trim(maskfile) // ' does not exist'
              stop
           else
              write (*,'(a)') '  mask == ' // trim(maskfile)
           end if
        else
           CALL getArgument(2, argument)
           mapfile2 = TRIM(adjustl(argument))
           INQUIRE(file=mapfile2, exist=there)
           IF (.NOT. there) STOP 'mapfile2 does not exist!'
           twomaps = .true.
           write (*,'(a)') ' Got two maps'
           WRITE (*,'(a)') ' Map2 == ' // TRIM(mapfile2)
        end if
        iarg = iarg + 1
     end do
  END IF

  if (twomaps) then
     call getmap(mapfile, map1, npix, nside, nstokes)
     call getmap(mapfile2, map2, npix2, nside2, nstokes2)
     nstokes = minval((/ nstokes, nstokes2 /))
     call fitmap(map1, map2, npix, nstokes)
     call diffmap(map, map1, map2, npix, nstokes)
  else
     call getmap(mapfile, map, npix, nside, nstokes)
  end if

  write (*,'("    npix == ",i0)') npix
  write (*,'(" nstokes == ",i0)') nstokes
  
  if (len_trim(maskfile) > 0) then
     call getmask(maskfile, mask, npix, nside)
  else
     allocate(mask(0:npix-1), stat=ierr)
     if (ierr /= 0) stop 'No room for mask'
     mask = .true. ! Start with all pixels included
  end if

  where(map(:,1) == hpx_dbadval) mask = .false.

  !do ipix = 0,npix-1
  !   if (abs(map(ipix,1)) < 1e-30) mask(ipix) = .false.
  !end do

  nhit = count(mask)

  write (*,'(/,"Sky coverage : ",i0," / ",i0," = ",f5.3)') &
       nhit, npix, dble(nhit)/npix

  do istokes = 1, nstokes
     if (nstokes == 1 .or. nstokes == 3) then
        write (*, '(/,"Stokes ",a," stats")') stokes(istokes)
     else
        write (*, '(/,"Column ",i0," stats")') istokes
     end if

     ! maxval and minval fail with large maps when using a mask
     mapmin = 1e30
     mapmax = -1e30
     do ipix = 0, npix-1
        if (.not. mask(ipix)) cycle
        if (map(ipix,istokes) > mapmax) mapmax = map(ipix,istokes)
        if (map(ipix,istokes) < mapmin) mapmin = map(ipix,istokes)
     end do

     write (*, '("Min   == ",es15.8,"  ")',advance='no') mapmin
     write (*, '("Max   == ",es15.8)') mapmax

     !mean = sum(map(:,istokes)) / nhit
     mean = sum(map(:, istokes), mask=mask) / nhit
     ! F90 where segfaults with nside2048 maps...
     !do ipix = 0,npix-1
     !   if (mask(ipix)) map(ipix,istokes) = map(ipix,istokes) - mean
     !end do
     
     stdev = sqrt(sum((map(:, istokes)-mean)**2, mask=mask) / nhit)
     !stdev = sqrt(dot_product(map(:,istokes), map(:,istokes)) / nhit)

     write (*, '("Mean  == ",es15.8,"  ")', advance='no') mean
     write (*, '("Stdev == ",es15.8)') stdev
  end do



contains



  subroutine getmap(mapfile, map, npix, nside, nstokes)
    character(len=*), intent(in) :: mapfile
    real(dp), pointer, intent(out) :: map(:,:)
    integer(i8b), intent(out) :: npix
    integer(i4b), intent(out) :: nside, nstokes

    integer(i4b) :: ordering

    ! Read in the map
    CALL tic
    npix = getsize_fits(mapfile, nmaps=nstokes, ordering=ordering, &
         nside=nside)
    if (associated(map)) deallocate(map)
    ALLOCATE(map(0:npix-1,nstokes), stat=ierr)
    IF (ierr /= 0) STOP 'No room for maps'
    CALL input_map(mapfile, map, npix, nstokes, fmissval=hpx_dbadval)
    IF (ordering == 1) THEN
       ! Covariance matrix is in NESTED scheme, so must the map be
       WRITE (*,*) ' Note: converting map into NESTED scheme'
       CALL convert_ring2nest(nside, map)
    END IF
    CALL toc('read map')

  end subroutine getmap



  subroutine getmask(maskfile, mask, npix, nside)
    character(len=*), intent(in) :: maskfile
    logical, pointer, intent(out) :: mask(:)
    integer(i8b), intent(in) :: npix
    integer(i4b), intent(in) :: nside

    real(dp), pointer :: masktemp(:,:)
    integer(i4b) :: ordering, nside_mask, nstokes_mask
    integer(i8b) :: npix_mask

    ! Read in the mask
    CALL tic
    npix_mask = getsize_fits(maskfile, nmaps=nstokes_mask, ordering=ordering, &
         nside=nside_mask)
    if (associated(mask)) deallocate(mask)
    ALLOCATE(masktemp(0:npix_mask-1,nstokes_mask), mask(0:npix-1), stat=ierr)
    IF (ierr /= 0) STOP 'No room for mask'
    CALL input_map(maskfile, masktemp, npix_mask, nstokes_mask)
    IF (ordering == 1) THEN
       ! Covariance matrix is in NESTED scheme, so must the map be
       WRITE (*,*) ' Note: converting map into NESTED scheme'
       CALL convert_ring2nest(nside_mask, masktemp)
    END IF
    CALL toc('read mask')

    if (nside /= nside_mask) then
       write (*,'(a,i0,a,i0)') 'ERROR: changing mask resolution not implemented: ', &
            nside, ' != ', nside_mask
       stop
    endif

    mask = (masktemp(:,0) > 1e-30)
    deallocate(masktemp)

    write (*,'(a,f5.2,a)') &
         'Mask excludes ',count(.not. mask) / dble(npix) * 100, ' % of the sky'

  end subroutine getmask



  subroutine diffmap(diff, map1, map2, npix, nstokes)
    real(dp), pointer, intent(out) :: diff(:,:)
    real(dp), pointer, intent(in) :: map1(:,:), map2(:,:)
    integer(i8b), intent(in) :: npix
    integer(i4b), intent(in) :: nstokes

    if (associated(diff)) deallocate(diff)
    ALLOCATE(diff(0:npix-1,nstokes), stat=ierr)
    IF (ierr /= 0) STOP 'No room for diff map'

    where (map1(:,1:nstokes) /= hpx_dbadval .and. map2(:,1:nstokes) /= hpx_dbadval)
       diff = map1(:,1:nstokes) - map2(:,1:nstokes)
    elsewhere
       diff = hpx_dbadval
    end where

  end subroutine diffmap



  subroutine fitmap(map1, map2, npix, nstokes)
    !
    ! Fit map1 = a + b * map2 using linear regression
    !
    real(dp), pointer, intent(in) :: map1(:,:), map2(:,:)
    integer(i8b), intent(in) :: npix
    integer(i4b), intent(in) :: nstokes

    integer(i8b) :: ipix
    integer(i4b) :: istokes
    real(dp) :: cov(2, 2), invcov(2, 2), dotp(2), coef(2)

    cov = 0

    do istokes = 1,nstokes
       do ipix = 0,npix-1
          if (map1(ipix, istokes) /= hpx_dbadval &
               .and. map2(ipix, istokes) /= hpx_dbadval) then
             cov(1, 1) = cov(1, 1) + map2(ipix, istokes)**2
             cov(2, 1) = cov(2, 1) + map2(ipix, istokes)
             cov(2, 2) = cov(2, 2) + 1
          endif
       end do
    end do
    cov(1, 2) = cov(2, 1)
    invcov(1, 1) = cov(2, 2)
    invcov(2, 1) = -cov(1, 2)
    invcov(1, 2) = invcov(2, 1)
    invcov(2, 2) = cov(1, 1)
    invcov = invcov / (cov(1, 1)*cov(2, 2) - cov(2, 1)*cov(1, 2))

    dotp = 0
    do istokes = 1,nstokes
       do ipix = 0,npix-1
          if (map1(ipix, istokes) /= hpx_dbadval .and. map2(ipix, istokes) /= hpx_dbadval) then
             dotp(1) = dotp(1) + map1(ipix, istokes) * map2(ipix, istokes)
             dotp(2) = dotp(2) + map1(ipix, istokes)
          endif
       end do
    end do

    coef(1) = dot_product(invcov(:, 1), dotp)
    coef(2) = dot_product(invcov(:, 2), dotp)
        
    write (*,'(a,g12.5,a,g12.5,a)') 'map1 = ',coef(2),' + ',coef(1),' * map2'

  end subroutine fitmap

END PROGRAM mapstats
