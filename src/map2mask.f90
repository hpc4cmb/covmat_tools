! Use pixel values and thresholds to construct a binary mask
!
! 07-16-2012 RK

PROGRAM map2mask

  USE healpix_types
  USE fitstools, ONLY   : output_map, getsize_fits, input_map, input_tod
  USE head_fits, ONLY   : add_card, write_minimal_header
  USE pix_tools, ONLY   : convert_ring2nest, nside2npix
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument
  CHARACTER(len=filenamelen) :: file_map, file_mask
  CHARACTER(len=80) :: header(1000)
  INTEGER(i4b) :: nmaps, ncols, nside, iargument, ordering, ierr
  integer(i8b) :: npix_tot, npix
  REAL(sp), allocatable :: map(:,:), mask(:,:)
  REAL(dp) :: maskmin, maskmax

  IF (nArguments() == 0)  THEN
     write (*,'(/,a,/)') ' Usage: map2mask <file_map> ' // &
          '--min <min value> --max <max_value> ' // &
          '[-o <file_mask>]'
     stop
  ELSE
     CALL getArgument(1, file_map)
     WRITE (*,*) ' file_map = '//trim(file_map)

     file_mask = '!mask.fits'
     iargument = 2
     maskmin = -1e30
     maskmax = 1e30
     do
        if (nArguments() < iargument) exit
        CALL getArgument(iargument, argument)

        if (index(argument, '-min') > 0) then
           iargument = iargument + 1
           call getArgument(iargument, argument)
           read (argument, *, iostat=ierr) maskmin
           if (ierr /= 0) then
              print *,'ERROR: could not parse min value from ' // trim(argument)
              stop
           end if
           WRITE (*,*) ' min value = ', maskmin
        else if (index(argument, '-max') > 0) then
           iargument = iargument + 1
           call getArgument(iargument, argument)
           read (argument, *, iostat=ierr) maskmax
           if (ierr /= 0) then
              print *,'ERROR: could not parse max value from ' // trim(argument)
              stop
           end if
           WRITE (*,*) ' max value = ', maskmax
        else if (index(argument, '-o') > 0) then
           iargument = iargument + 1
           call getArgument(iargument, file_mask)
        else
           write (*,*) 'ERROR: unrecognized option: '//trim(argument)
           stop
        end if

        iargument = iargument + 1
     end do

     WRITE (*,*) ' Generating ' // TRIM(file_mask)
  END IF

  npix_tot = getsize_fits(file_map, nmaps=ncols, ordering=ordering, nside=nside)
  if (ncols > 1) print *,'Warning, only examining the first column out of ',ncols
  nmaps = ncols
  npix = nside2npix(nside)

  write (*,'(" nmaps == ",i8)') nmaps
  write (*,'(" nside == ",i8)') nside
  write (*,'("  npix == ",i8)') npix

  ALLOCATE(map(0:npix-1,1), mask(0:npix-1,1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for map or mask'

  CALL tic
  call input_map(file_map, map, npix, 1, fmissval=HPX_SBADVAL)
  call toc('Read map')
  
  call tic
  where (map >= maskmin .and. map <= maskmax)
     mask = 1
  elsewhere
     mask = 0
  end where
  write (*,'(a,f6.3)') ' Sky fraction = ',dble(sum(mask))/npix
  CALL toc('Generate mask')

  CALL write_minimal_header(header, 'MAP', nside=nside, order=ordering, &
       creator='map2mask', polar=.false.)

  CALL output_map(mask, header, file_mask)


END PROGRAM map2mask
