! This code expands a sparse covariance matrix into a full
! matrix according to the pixel list file


PROGRAM sparse2full

  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, file_covmat, outfile, pixellistfile
  CHARACTER(len=80)          :: spixel
  INTEGER(I8B)               :: pixel,npix,nstokes,irow,nside,col,istokes
  INTEGER(I8B)               :: colOut
  INTEGER, PARAMETER         :: covmat_unit=55
  INTEGER(dp), PARAMETER     :: elements_max=1E8
  REAL(dp), PARAMETER        :: not_read=-1E30
  REAL(dp)                   :: row(elements_max)
  REAL(dp), POINTER          :: covmat(:,:), temp(:)
  INTEGER, POINTER           :: hit_pixels(:), pixellist(:)
  LOGICAL                    :: there, reduce=.TRUE.
  LOGICAL                    :: logscale=.FALSE., normalize=.FALSE.
  LOGICAL                    :: interlaced=.FALSE., use_pixellist=.FALSE.
  !LOGICAL                    :: logscale=.TRUE., normalize=.TRUE.
  INTEGER(I4B)               :: ierr, isig, iSig2, iArgument

  IF (nArguments() < 4)  THEN
     WRITE (*,*) 'Usage: sparse2full <covmat> <pixellist> ' // &
          '<nstokes> [--interlaced] [--out outfile] '
     STOP
  ELSE
     CALL getArgument(1, argument)
     file_covmat = TRIM(argument)
     outroot     = file_covmat
     INQUIRE(file=TRIM(file_covmat), exist=there)
     IF (.NOT. there) STOP 'file_covmat does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(file_covmat)

     CALL getArgument(2, argument)
     INQUIRE(file=TRIM(argument), exist=there)
     if (there) then
        ! user provided pixel list
        pixellistfile = trim(argument)
        write (*,'(a)') ' Mapping pixels according to '//trim(pixellistfile)
        use_pixellist = .true.
     else
        stop 'pixellist file does not exist!'
     end if

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i8)') 'Nstokes == ', nstokes

     iArgument = 4
     outfile = 'full_'//file_covmat
     DO
        IF (iArgument > nArguments()) EXIT

        CALL getArgument(iArgument, argument)
        IF (INDEX(argument, '--interlaced') /= 0) THEN
           WRITE (*,'(a)') ' Assuming input is in interlaced form'
           interlaced = .TRUE.
        ELSE IF (INDEX(argument, '-o') /= 0) THEN
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           outfile = TRIM(argument)
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

  IF (INDEX(file_covmat, '.fits') /= 0) THEN
     ! fits file
     IF (pixel < 0) STOP 'pixel < 0 not implemented for fits files'
     DO istokes = 1, nstokes
        col  = MODULO(pixel, npix) + npix*(istokes-1)
        temp = get_col_from_fits_matrix(col, npix, nstokes, file_covmat)
        covmat(:, istokes) = temp
     END DO
  ELSE
     ! File is in C binary format
     OPEN(unit=covmat_unit, file=TRIM(file_covmat), status='old', &
          form='unformatted', access='direct', recl=npix*nstokes*8)
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


  DEALLOCATE(covmat)

END PROGRAM sparse2full
