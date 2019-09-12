! This code removes from the covariance matrix rows
! and colums excluded by the supplied mask

PROGRAM mask_covmat

  USE healpix_types
  use udgrade_nr, only  : udgrade_nest
  USE fitstools, ONLY   : input_map, getsize_fits
  use pix_tools, only   : convert_ring2nest
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, mask_file, out_file
  CHARACTER(len=filenamelen) :: mask2_file
  INTEGER(i4b) :: nstokes, nside
  INTEGER(i8b) :: icol, icol_out, npix_temp, npix_pol, npix
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:)
  real(dp), pointer :: mask_temp(:), mask_pol(:), mask(:,:)
  logical, pointer :: lmask(:)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, ordering, iArgument, rec_len

  IF (nArguments() < 4)  THEN
     WRITE (*,'(/,a,/)') '  Usage: mask_covmat <covmat> <nside> <nmaps> ' // &
          '<mask> [<mask2>] [-o <outfile>]'
     STOP
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

     CALL getArgument(4, argument)
     mask_file = trim(argument)
     INQUIRE(file=TRIM(mask_file), exist=there)
     IF (.NOT. there) STOP 'mask does not exist'
     WRITE (*,*) ' Reading ' // TRIM(mask_file)

     out_file   = 'masked_'//TRIM(ADJUSTL(covmat_file))
     mask2_file = ''
     iArgument = 5
     do
        if (nArguments() < iArgument) exit

        CALL getArgument(iArgument, argument)
        if (index(argument,'-o') /= 0) then
           CALL getArgument(iArgument+1, argument)
           out_file = TRIM(ADJUSTL(argument))
           iArgument = iArgument + 1
        else
           mask2_file = trim(argument)
           INQUIRE(file=TRIM(mask2_file), exist=there)
           IF (.NOT. there) STOP 'mask2 does not exist'
           write (*,*) 'Reading separate mask for polarization: ' &
                // TRIM(mask2_file)
        end if

        iArgument = iArgument + 1
     end do

     WRITE (*,*) ' Writing to ' // TRIM(out_file)
  END IF

  npix = 12*nside**2
  allocate(mask(npix,2), stat=ierr)
  if (ierr /= 0) stop 'No room for masks'

  ! Read in the mask(s)
  CALL tic
  call read_mask(mask(:,1), mask_file, nside)
  CALL toc('read mask')

  ! read the possible separate polarization mask
  if (trim(mask2_file) == '') then
     mask(:,2) = mask(:,1)
  else
     call tic
     call read_mask(mask(:,2), mask2_file, nside)
     CALL toc('read mask2')
  end if

  ! final parsing
  mask_temp => mask(:,1)
  mask_pol  => mask(:,2)
  npix_temp = count(mask_temp /= 0)
  npix_pol  = count(mask_pol /= 0)
  write (*,'(a,i0,"+",i0,a)') '  Mask leaves ', npix_temp, npix_pol, ' pixels'


  if (nstokes == 1) npix_pol = 0
  ALLOCATE(covmat(0:npix*nstokes-1), lmask(0:npix*nstokes+1), stat=ierr)

  ! finally put together a mask for covmat columns
  if (nstokes == 1) then
     lmask    = mask_temp /= 0
     npix_pol = 0
  else
     lmask(     0:  npix-1) = mask_temp /= 0
     lmask(  npix:2*npix-1) = mask_pol  /= 0
     lmask(2*npix:3*npix-1) = mask_pol  /= 0
  end if

  IF (ierr /= 0) STOP 'No room for covmat'
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  inquire(iolength=rec_len) pack(covmat, lmask)
  OPEN(unit=covmat_unit+1, file=TRIM(out_file), status='replace', &
       form='unformatted', access='direct', recl=rec_len)

  ! read, mask and write. one column at a time
  CALL tic
  icol_out = 0
  DO icol = 0, npix*nstokes-1
     if (lmask(icol)) then
        READ(covmat_unit, rec=icol+1) covmat
        icol_out = icol_out + 1
        WRITE(covmat_unit+1, rec=icol_out) pack(covmat, lmask)
     end if
  END DO
  CLOSE(covmat_unit)
  CLOSE(covmat_unit+1)
  call toc('Mask covmat')

  ! Write an ASCII file that describes mapping between
  ! pixel indices and pixel numbers
  covmat = (/ (modulo(icol,npix), icol=0, npix*nstokes-1) /)
  open(unit=covmat_unit+2, file='index2pixel.dat', status='replace', &
       form='formatted')
  write (covmat_unit+2,*) npix_temp, npix_pol
  write (covmat_unit+2,'(10000000(i0,/))') int(pack(covmat, lmask), i4b)
  write (*,*) ' Pixel indices saved in index2pixel.dat'

  CLOSE(covmat_unit+2)
  CALL toc('mask_covmat')


contains


  subroutine read_mask(mask, mask_file, nside)

    real(dp) :: mask(:)
    character(len=filenamelen) :: mask_file
    integer(i4b) :: nside

    integer(i4b) :: nmask, nside_mask
    integer(i8b) :: npix_mask
    real(dp), pointer :: mask_in(:,:)

    npix = getsize_fits(mask_file, nmaps=nmask, ordering=ordering, &
         nside=nside_mask)
    if (nmask > 1) write (*,'(a,a,i0,a)' )'WARNING: ', trim(mask_file), &
         ' contains ', nmask, ' columns, using only first'
    npix_mask = npix

    ALLOCATE(mask_in(0:npix_mask-1, nmask), stat=ierr)
    IF (ierr /= 0) STOP 'No room for mask'

    CALL input_map(mask_file, mask_in, npix_mask, nmask)
    IF (ordering == 1) THEN
       ! Covariance matrix is in NESTED scheme, so must the mask be
       WRITE (*,'(a,a,a)') ' Note: converting ', trim(mask_file), &
            ' into NESTED scheme'
       CALL convert_ring2nest(nside_mask, mask_in)
    END IF
    
    ! possibly downgrade the mask
    if (nside == nside_mask) then
       mask = mask_in(:,1)
    else
       if (ierr /= 0) stop 'No room for mask'
       call udgrade_nest(mask_in(:,1), nside_mask, mask, nside)
    end if

    deallocate(mask_in)
  end subroutine read_mask

END PROGRAM mask_covmat
