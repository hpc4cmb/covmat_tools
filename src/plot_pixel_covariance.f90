! This code reads madam pixel-pixel covariance matrix
! and plots the given pixel covariance maps

PROGRAM plot_pixel_covariance

  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc, get_matrix_size_dp

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outroot
  CHARACTER(len=80) :: spixel
  INTEGER(I8B) :: nelem, npix, pixel, col
  integer(i4b) :: nstokes, istokes, nside
  INTEGER, PARAMETER :: covmat_unit=55
  INTEGER(i8b), PARAMETER :: elements_max=1E8
  REAL(dp), PARAMETER :: not_read=-1E30
  REAL(dp), POINTER :: covmat(:, :), temp(:)
  LOGICAL :: there
  LOGICAL :: logscale=.FALSE., normalize=.FALSE.
  LOGICAL :: mask_diagonal=.FALSE.
  LOGICAL :: interlaced=.FALSE.
  INTEGER(I4B) :: ierr, isig, iSig2, iArgument, rec_len

  IF (nArguments() < 2)  THEN
     WRITE (*,*) 'Usage: plot_pixel_covariance <covmat>' // &
          ' <row> [--interlaced] [--normalize]' // &
          ' [--mask_diagonal] [--out <outroot>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     outroot     = covmat_file
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, spixel)
     READ (spixel, *, iostat=ierr) pixel
     if (ierr /= 0) stop 'Unable to parse pixel'
     IF (pixel < 0) THEN
        WRITE (*,*) 'Row is negative! Plotting the diagonal instead.'
     ELSE
        WRITE (*,'(a,i8)') 'pixel   == ', pixel
     END IF

     iArgument = 3
     DO
        IF (iArgument > nArguments()) EXIT

        CALL getArgument(iArgument, argument)
        IF (INDEX(argument, '--interlaced') /= 0) THEN
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

  call get_matrix_size_dp(covmat_file, nside, nstokes, npix, nelem)
  WRITE (*,'(a,i8)') 'Nside   == ', nside
  WRITE (*,'(a,i8)') 'Nstokes == ', nstokes
  npix = npix / nstokes

  ALLOCATE(covmat(0:npix*nstokes-1,nstokes), temp(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat'

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
        if (interlaced) then
           !stop 'pixel < 0 not implemented for interlaced files'
           !DO col = 0, npix-1
           !   DO isig = 1, nstokes
           !      READ(covmat_unit, rec=isig+col*nstokes) temp
           !      covmat(col:npix*nstokes-1:npix,isig) = &
           !           temp(modulo(col,npix)+isig-1:npix*nstokes-1:npix)
           !   END DO
           !end DO
           DO col = 0, npix-1
              DO isig = 1, nstokes
                 READ(covmat_unit, rec=col*nstokes+isig) temp
                 covmat(col+(isig-1)*npix,:) = temp(col*nstokes:(col+1)*nstokes-1)
              end DO
           END DO
        else
           DO col = 0, npix*nstokes-1
              READ(covmat_unit, rec=col+1) temp
              covmat(col,:) = temp(modulo(col,npix):npix*nstokes-1:npix)
           END DO
        end if
     ELSE
        ! Plot the columns corresponding to requested pixel
        DO isig = 0, nstokes-1
           IF (interlaced) THEN
              READ(covmat_unit, rec=pixel*nstokes+isig+1) temp
              DO isig2 = 0, nstokes-1
                 covmat(npix*isig2:npix*(isig2+1)-1, isig+1) = &
                      temp(isig2:nstokes*npix-1:nstokes)
              END DO
           ELSE
              READ(covmat_unit, rec=pixel+isig*npix+1) temp
              covmat(:, isig+1) = temp
           END IF
        END DO
        IF (normalize) then
           write (*,'(a)', advance='no') 'Normalizing by pixel stdevs.. '
           ! normalize by reference pixel stdev
           do isig = 1, nstokes
              if (covmat(pixel+(isig-1)*npix, isig).gt.0) then
                covmat(:, isig) = covmat(:, isig) / &
                   sqrt(covmat(pixel+(isig-1)*npix, isig))
              else
                 covmat(:, isig) = HPX_DBADVAL
              endif
           end do
           ! normalize each pixel by corresponding stdev (sqrt(diag))
           DO col = 0, npix*nstokes-1
              READ(covmat_unit, rec=col+1) temp
              do isig = 1, nstokes
                 if (temp(col).gt.0.and.covmat(col,isig).gt.HPX_DBADVAL) then
                    covmat(col,isig) = covmat(col,isig) / sqrt(temp(col))
                 else
                    covmat(col,isig) = HPX_DBADVAL
                 endif
              enddo
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

    integer(i8b) :: icol, npix
    integer(i4b) :: nstokes
    CHARACTER(filenamelen) :: covmat_file

    REAL(dp) :: get_col_from_fits_matrix(0:npix*nstokes-1)

    REAL(dp), POINTER :: temp(:, :)
    INTEGER(I4B) :: ierr

    ALLOCATE(temp(0:npix*nstokes-1, 1), stat=ierr)
    IF (ierr /= 0) STOP 'Out of memory in get_col_from_fits_matrix'

    CALL input_tod(covmat_file, temp, int(npix*nstokes, i8b), 1, &
         firstpix=int(icol*npix*nstokes, i8b))

    get_col_from_fits_matrix = temp(:,1)

  END FUNCTION get_col_from_fits_matrix



  SUBROUTINE write_single_covariance_map(covmat, pixel, outfileroot)
    ! Writes the covariance of the given pixel into file
    ! the array covmat is dimension(nstokes*npix, nstokes) and
    ! contains the three columns of the covariance matrix that
    ! correspond to the reference pixel in nstokes maps

    USE fitstools, ONLY : output_map
    USE head_fits, ONLY : add_card, write_minimal_header

    IMPLICIT NONE

    REAL(dp), POINTER :: covmat(:, :)
    INTEGER(i8b) :: pixel
    CHARACTER(len=*) :: outfileroot

    CHARACTER(len=filenamelen) :: outfitsfile
    INTEGER :: sig1, sig2
    INTEGER(i8b) :: pix
    REAL(dp), POINTER :: map(:, :)
    CHARACTER(len=80), DIMENSION(1000) :: header

    REAL(dp), POINTER :: row(:)

    ALLOCATE(row(0:npix*nstokes-1))
    ALLOCATE(map(0:npix-1, 1:nstokes))

    pixel = MODULO(pixel, npix)

    DO sig1 = 1,nstokes

       WRITE (*,*) 'Processing row #', pixel

       row   = covmat(:,sig1)

       !write (*,*) row(pixel:pixel+5)
       IF (COUNT(row/=0)==0) WRITE (*,'(a,i10,a)') &
            'Warning! map for pixel # ', pixel, ' is empty!'

       IF (logscale) THEN
          WRITE (*,*) 'Writing map in logscale'
          DO pix = 0, npix*nstokes-1
             IF (row(pix)/=0) row(pix) = LOG(ABS(row(pix)))
          END DO
       END IF

       WRITE (spixel,'(i10)') pixel

       outfitsfile='!'//TRIM(outfileroot)//'_row_'//&
            TRIM(ADJUSTL(spixel))//'.fits'
       !outfitsfile = '!row_'//TRIM(ADJUSTL(spixel))//'.fits'
       WRITE (*,*) 'Writing to : ' // TRIM(outfitsfile) // ' ...'

       DO sig2 = 1, nstokes
          map(:,sig2) = row(npix*(sig2-1):npix*sig2-1)
       END DO

       !IF (normalize) THEN
       !   WRITE (*,*) 'Normalizing diagonal to one'
       !   norm = map(MODULO(pixel,npix), sig1)
       !   map  = map/norm
       !END IF

       IF (mask_diagonal) THEN
          ! set actual pixel variances to zero to bring out more structure
          map (MODULO(pixel,npix), :) = 0.0
       END IF       

       header = ''
       call write_minimal_header(header, 'MAP', nside=int(nside,i4b), &
            ordering='Nested', creator='plot_pixel_covariance', polar=nstokes==3)            

       CALL output_map(map, header, outfitsfile)

       pixel = pixel + npix

    END DO

  END SUBROUTINE write_single_covariance_map



  SUBROUTINE write_single_variance_map(covmat, outfileroot)
    ! Writes the diagonal of the covariance matrix as a map

    USE fitstools, ONLY : output_map
    USE head_fits, ONLY : add_card, write_minimal_header

    IMPLICIT NONE

    REAL(dp), POINTER :: covmat(:,:)
    CHARACTER(len=*) :: outfileroot

    CHARACTER(len=filenamelen) :: outfitsfile
    INTEGER :: sig1, sig2
    INTEGER(dp) :: pix
    REAL(dp), POINTER :: map(:,:)
    CHARACTER(len=80), DIMENSION(1000) :: header

    ALLOCATE(map(0:npix-1, 1:nstokes))

    header = ''
    call write_minimal_header(header, 'MAP', nside=int(nside,i4b), &
            ordering='Nested', creator='plot_pixel_covariance', polar=nstokes==3)            

    do sig2 = 1, nstokes
       DO sig1 = 1, nstokes
          map(:,sig1) = covmat((sig1-1)*npix:sig1*npix-1, sig2)
          IF (logscale) THEN
             WRITE (*,*) 'Writing map in logscale'
             DO pix = 0, npix-1
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



END PROGRAM plot_pixel_covariance
