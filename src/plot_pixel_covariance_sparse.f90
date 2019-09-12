! This code reads madam pixel-pixel covariance matrix
! and plots the given pixel covriance maps

PROGRAM plot_pixel_covariance_sparse

  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outroot, pixellistfile
  CHARACTER(len=80) :: spixel
  INTEGER(I4B) :: nstokes, nside, istokes
  INTEGER(I8B) :: colOut, pixel, npix, col
  INTEGER, PARAMETER :: covmat_unit=55
  INTEGER(dp), PARAMETER :: elements_max=1E8
  REAL(dp), PARAMETER :: not_read=-1E30
  REAL(dp), POINTER :: covmat(:,:), temp(:)
  INTEGER, POINTER :: pixellist(:)
  LOGICAL :: there
  LOGICAL :: logscale=.FALSE., normalize=.FALSE.
  LOGICAL :: mask_diagonal=.FALSE., cMatrix=.TRUE.
  LOGICAL :: interlaced=.FALSE., use_pixellist=.FALSE.
  INTEGER(I4B) :: ierr, isig, iSig2, iArgument, rec_len

  IF (nArguments() < 4)  THEN
     WRITE (*,*) 'Usage: plot_pixel_covariance <covmat> <pixellist> ' // &
          '<nstokes> <row> [--F] [--interlaced] [--normalize] ' // &
          '[--mask_diagonal] [--out <outroot>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     outroot     = covmat_file
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     INQUIRE(file=TRIM(argument), exist=there)
     if (there) then
        ! user provided pixel list
        pixellistfile = trim(argument)
        write (*,'(a)') ' Mapping pixels according to '//trim(pixellistfile)
        use_pixellist = .true.
     else
        READ (argument, *) nside
     WRITE (*,'(a,i8)') 'Nside   == ', nside
     end if

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i8)') 'Nstokes == ', nstokes

     CALL getArgument(4, spixel)
     READ (spixel, *) pixel
     IF (pixel < 0) THEN
        WRITE (*,*) 'Row is negative! Plotting the diagonal instead.'
     ELSE
        WRITE (*,'(a,i8)') 'pixel   == ', pixel
     END IF

     iArgument = 5
     DO
        IF (iArgument > nArguments()) EXIT

        CALL getArgument(iArgument, argument)
        IF (INDEX(argument, '-F') /= 0) THEN
           WRITE (*,'(a)') ' Assuming input is Fortran binary'
           cMatrix = .FALSE.
        ELSE IF (INDEX(argument, '--interlaced') /= 0) THEN
           WRITE (*,'(a)') ' Assuming input is in interlaced form'
           interlaced = .TRUE.
        ELSE IF (INDEX(argument, '--normalize') /= 0) THEN
           WRITE (*,'(a)') ' Normalizing row with diagonal amplitude'
           normalize = .TRUE.
        ELSE IF (INDEX(argument, '--mask_diagonal') /= 0) THEN
           WRITE (*,'(a)') ' Masking the diagonal elements by zero'
           mask_diagonal = .TRUE.
        ELSE IF (INDEX(argument, '-o') /= 0) THEN
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           outroot = TRIM(argument)
        ELSE
           WRITE (*,*) ' Unrecognized argument: ' // TRIM(ADJUSTL(argument))
           STOP
        END IF
        iArgument = iArgument + 1
     END DO
  END IF

  WRITE (*,'(a)') ' Using '//TRIM(outroot)//' as root name for output'

  if (use_pixellist) then
     ! first line of a pixellistfile lists the nside and number of hit pixels
     open(file=pixellistfile,unit=55,status='old')
     read(55,*) nside, npix
     write (*,'(a,i0)') ' nside == ', nside
     write (*,'(a,i0)') ' npix == ', npix
     allocate(pixellist(0:npix*nstokes-1), stat=ierr)
     if (ierr /= 0) stop 'No room for pixellist'
     read(55,*) pixellist(0:npix-1)
     close(55)
     if (nstokes == 3) then
        pixellist(npix:2*npix-1) = pixellist(0:npix-1)+12*nside**2
        pixellist(2*npix:3*npix-1) = pixellist(0:npix-1)+2*12*nside**2
     end if
  else
     npix = 12*nside**2
  end if

  ! output has 12*nside**2 pixels, but input (temp) must conform to pixellist
  ALLOCATE(covmat(0:12*nside**2*nstokes-1,nstokes), temp(0:npix*nstokes-1), &
       stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat'
  covmat = hpx_dbadval ! initialize to null value

  CALL tic

  IF (INDEX(covmat_file, '.fits') /= 0) THEN
     ! fits file
     IF (pixel < 0) STOP 'pixel < 0 not implemented for fits files'
     DO istokes = 1, nstokes
        col  = MODULO(pixel, npix) + npix*(istokes-1)
        temp = get_col_from_fits_matrix(col, npix, nstokes, covmat_file)
        covmat(:, istokes) = temp
     END DO
  ELSE
     ! File is in C binary format
     inquire(iolength=rec_len) temp
     OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
          form='unformatted', access='direct', recl=rec_len)
     IF (pixel < 0) THEN
        ! plot the diagonal
        DO col = 0, npix*nstokes-1
           READ(covmat_unit, rec=col+1) temp
           colOut = col
           if (use_pixellist) colOut = pixellist(col)
           covmat(colOut,:) = temp(modulo(col,npix):npix*nstokes-1:npix)
        END DO
     ELSE
        ! Plot the columns corresponding to requested pixel
        DO isig = 0, nstokes-1
           IF (interlaced) THEN
              READ(covmat_unit, rec=pixel*nstokes+isig+1) temp
              DO isig2 = 0, nstokes-1
                 if (use_pixellist) then
                    covmat(pixellist(npix*isig2:npix*(isig2+1)-1), isig+1) = &
                         temp(isig2:nstokes*npix-1:nstokes)
                 else
                    covmat(npix*isig2:npix*(isig2+1)-1, isig+1) = &
                         temp(isig2:nstokes*npix-1:nstokes)
                 end if
              END DO
           ELSE
              READ(covmat_unit, rec=pixel+isig*npix+1) temp
              if (use_pixellist) then
                 covmat(pixellist, isig+1) = temp
              else
                 covmat(:, isig+1) = temp
              end if
           END IF
        END DO
        IF (normalize) then
           write (*,'(a)', advance='no') 'Normalizig by pixel stdevs.. '
           ! normalize by reference pixel stdev
           do isig = 1, nstokes
              covmat(:, isig) = covmat(:, isig) / &
                   sqrt(covmat(pixel+(isig-1)*npix, isig))
           end do
           ! normalize each pixel by corresponding stdev (sqrt(diag))
           DO col = 0, npix*nstokes-1
              READ(covmat_unit, rec=col+1) temp
              covmat(col,:) = covmat(col,:) / sqrt(temp(col))
           END DO
           write (*,*) 'done!'
        end IF
     END IF
     CLOSE(unit=covmat_unit)
  END IF

  CALL toc('Read matrix')

  IF (pixel < 0) THEN
     CALL write_single_variance_map(covmat, covmat_file)
  ELSE
     if (use_pixellist) pixel = pixellist(pixel)
     CALL write_single_covariance_map(covmat, pixel, outroot)
  END IF

  DEALLOCATE(covmat)



CONTAINS


  FUNCTION get_col_from_fits_matrix(icol, npix, nstokes, covmat_file)
    ! This function returns the requested column of a covariance matrix
    ! from a double precision, single column fits file
    !
    ! February 27th 2008 Reijo Keskitalo

    USE healpix_types
    USE fitstools, ONLY : input_tod

    IMPLICIT NONE

    INTEGER(I8B) :: icol, npix
    INTEGER(I4B) :: nstokes
    CHARACTER(filenamelen) :: covmat_file

    REAL(dp) :: get_col_from_fits_matrix(0:npix*nstokes-1)

    REAL(dp), POINTER :: temp(:, :)
    INTEGER(I4B) :: ierr

    ALLOCATE(temp(0:npix*nstokes-1, 1), stat=ierr)
    IF (ierr /= 0) STOP 'Out of memory in get_col_from_fits_matrix'

    CALL input_tod(covmat_file, temp, npix*nstokes, 1, &
         firstpix=icol*npix*nstokes)

    get_col_from_fits_matrix = temp(:,1)

  END FUNCTION get_col_from_fits_matrix



  SUBROUTINE write_single_covariance_map(covmat, pixel, outfileroot)
    ! Writes the covariance of the given pixel into file
    ! the array covmat is dimension(nstokes*npix, nstokes) and
    ! contains the three columns of the covariance matrix that
    ! correspond to the reference pixel in nstokes maps

    USE fitstools, ONLY : output_map
    USE head_fits, ONLY : add_card

    IMPLICIT NONE

    REAL(dp), POINTER :: covmat(:, :)
    INTEGER(dp) :: pixel
    CHARACTER(len=*) :: outfileroot

    CHARACTER(len=1024) :: outfitsfile
    INTEGER :: sig1, sig2
    INTEGER(i8b) :: pix
    REAL(dp), POINTER :: map(:, :)
    CHARACTER(len=80), DIMENSION(1:10) :: header
    INTEGER :: itemp

    REAL(dp), POINTER :: row(:)

    ALLOCATE(row(0:12*nside**2*nstokes-1))
    ALLOCATE(map(0:12*nside**2-1, 1:nstokes))

    pixel = MODULO(pixel, int(12*nside**2, i8b))

    DO sig1 = 1,nstokes

       WRITE (*,*) 'Processing row #', pixel

       row = covmat(:, sig1)

       !write (*,*) row(pixel:pixel+5)
       IF (COUNT(row/=0) == 0) WRITE (*,'(a,i10,a)') &
            'Warning! map for pixel # ', pixel, ' is empty!'

       IF (logscale) THEN
          WRITE (*,*) 'Writing map in logscale'
          DO pix = 0, 12*nside**2*nstokes-1
             IF (row(pix)/=0) row(pix) = LOG(ABS(row(pix)))
          END DO
       END IF

       WRITE (spixel,'(i10)') pixel

       outfitsfile='!'//TRIM(outfileroot)//'_row_'//&
            TRIM(ADJUSTL(spixel))//'.fits'
       !outfitsfile = '!row_'//TRIM(ADJUSTL(spixel))//'.fits'
       WRITE (*,*) 'Writing to : ' // TRIM(outfitsfile) // ' ...'

       DO sig2 = 1, nstokes
          map(:,sig2) = row(12*nside**2*(sig2-1):12*nside**2*sig2-1)
       END DO

       !IF (normalize) THEN
       !   WRITE (*,*) 'Normalizing diagonal to one'
       !   norm = map(MODULO(pixel,npix), sig1)
       !   map  = map/norm
       !END IF

       IF (mask_diagonal) THEN
          ! set actual pixel variances to zero to bring out more structure
          map(MODULO(pixel, int(12*nside**2, i8b)), :) = 0
       END IF

       header(:) = ''
       CALL add_card(header, 'PIXTYPE',  'HEALPIX',  'HEALPIX Pixelisation')
       CALL add_card(header, 'ORDERING', 'NESTED',   'Pixel ordering scheme')
       itemp =  nside
       CALL add_card(header, 'NSIDE',     itemp,     'nside of the map')
       CALL add_card(header, 'FIRSTPIX',  0,         'First pixel # (0 based)')
       itemp = 12*nside**2-1
       CALL add_card(header, 'LASTPIX',   itemp,     'Last pixel # (0 based)')
       CALL add_card(header, 'POLAR',     nstokes==3,'Polarization flag')

       CALL output_map(map, header, outfitsfile)

       pixel = pixel + 12*nside**2

    END DO

  END SUBROUTINE write_single_covariance_map



  SUBROUTINE write_single_variance_map(covmat, outfileroot)
    ! Writes the diagonal of the covariance matrix as a map

    USE fitstools, ONLY : output_map
    USE head_fits, ONLY : add_card

    IMPLICIT NONE

    REAL(dp), POINTER :: covmat(:, :)
    CHARACTER(len=*) :: outfileroot

    CHARACTER(len=filenamelen) :: outfitsfile
    INTEGER(i4b) :: sig1, sig2
    INTEGER(i8b) :: pix
    REAL(dp), POINTER :: map(:, :)
    CHARACTER(len=80), DIMENSION(1:10) :: header
    INTEGER :: itemp

    ALLOCATE(map(0:12*nside**2-1, 1:nstokes))

    header(:) = ''
    CALL add_card(header, 'PIXTYPE',  'HEALPIX', 'HEALPIX Pixelisation')
    CALL add_card(header, 'ORDERING', 'NESTED',  'Pixel ordering scheme')
    itemp =  nside
    CALL add_card(header, 'NSIDE',     itemp,    'nside of the map')
    CALL add_card(header, 'FIRSTPIX',  0,        'First pixel # (0 based)')
    itemp = 12*nside**2-1
    CALL add_card(header, 'LASTPIX',   itemp,    'Last pixel # (0 based)')

    do sig2 = 1, nstokes
       DO sig1 = 1, nstokes
          map(:,sig1) = covmat((sig1-1)*12*nside**2:sig1*12*nside**2-1, sig2)
          IF (logscale) THEN
             WRITE (*,*) 'Writing map in logscale'
             DO pix = 0, 12*nside**2-1
                IF (map(pix,sig1)/=0) map(pix,sig1) = LOG(map(pix,sig1))
             END DO
          END IF
       END DO

       select case(sig2)
       case (1)
          ! Diagonals from II, IQ and IU blocks
          outfitsfile='!'//TRIM(outfileroot)//'_I_diagonals.fits'
       case (2)
          ! Diagonals from QI, QQ and QU blocks
          outfitsfile='!'//TRIM(outfileroot)//'_Q_diagonals.fits'
       case (3)
          ! Diagonals from UI, UQ and UU blocks
          outfitsfile='!'//TRIM(outfileroot)//'_U_diagonals.fits'
       end select
       WRITE (*,*) 'Writing to : ' // TRIM(outfitsfile) // ' ...'

       CALL output_map(map, header, outfitsfile)
    end do

  END SUBROUTINE write_single_variance_map



END PROGRAM plot_pixel_covariance_sparse
