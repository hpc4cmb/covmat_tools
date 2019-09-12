PROGRAM downgrade_evecs

  USE fitstools
  USE head_fits
  USE pix_tools
  USE healpix_types
  USE extension
  use covmat_util, only : get_matrix_size_dp

  use downgrade

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: file_in, file_out, file_nobs, file_beam
  INTEGER(i4b) :: nside_in, nside_out, midnside=-1
  INTEGER(i4b) :: ierr, imap, nmaps
  INTEGER(i8b) :: npix_in, npix_out, npixtot
  REAL(dp) :: FWHM, rcondlim
  CHARACTER(len=filenamelen) :: argument
  LOGICAL :: use_cosine_window=.FALSE., noiseweight=.false., use_beam=.false.
  logical :: smooth=.FALSE., polsmooth=.TRUE.

  CALL parseArguments()

  write (*,'(a)') 'Parsing size of ' // trim(file_in)
  call get_matrix_size_dp(file_in, nside_in, nmaps, npix_in, npixtot)
  WRITE (*,'(a,i8)') 'Nside_in == ', nside_in
  WRITE (*,'(a,i8)') 'Nmaps_in == ', nmaps

  CALL processEvecs()


CONTAINS


  SUBROUTINE processEvecs()

    ! Read, process and write the eigenvectors one by one

    REAL(dp), POINTER :: map_in(:, :), map_out(:, :), evec_in(:), evec_out(:), &
         map_mid(:, :), map_unsmoothed(:, :)
    INTEGER(i4b) :: ierr, evecUnit=55, rec_len
    INTEGER(i8b) :: icol, npix_mid
    
    if ( smooth ) then
       if (midnside > 0) then
          if (midnside > nside_in) stop 'midnside must be less than input nside'
          if (midnside < nside_out) stop 'midnside must be greater than or equal the output nside'
       else
          IF (nside_out < nside_in ) THEN
             midnside = nside_out * 2 
          ELSE
             midnside = nside_out
          END IF
       end if
    else
       if (midnside > 0) print *,'Warning, midnside ignored without smoothing'
       midnside = nside_out
    end if

    npix_in = nside2npix(nside_in)
    npix_mid = nside2npix(midnside)
    npix_out = nside2npix(nside_out)

    ALLOCATE(map_in(0:npix_in-1, nmaps), map_out(0:npix_out-1, nmaps), &
         map_mid(0:npix_mid-1, nmaps), map_unsmoothed(0:npix_out-1,nmaps), &
         evec_in(0:npix_in*nmaps-1), evec_out(0:npix_out*nmaps-1), stat=ierr)
    IF (ierr /= 0) STOP 'no room for evecs'

    inquire(iolength=rec_len) evec_in
    open(unit=evecUnit, file=trim(file_in), status='old', &
         form='unformatted', access='direct', recl=rec_len)

    inquire(iolength=rec_len) evec_out
    open(unit=evecUnit+1, file=trim(file_out), status='replace', &
         form='unformatted', access='direct', recl=rec_len)

    DO icol = 1, npix_in*nmaps
       IF (MODULO(icol, npix_in*nmaps/10) == 1) &
            WRITE (*,'(i4,a)') icol*100/(npix_in*nmaps), '% done'

       READ(evecunit, rec=icol) evec_in
       map_in = reshape(evec_in, (/npix_in, int(nmaps, i8b)/))

       if (nside_in /= midnside) then
          if (noiseweight) then
             call noise_weight_downgrade(file_nobs, map_in, map_mid, nside_in, &
                  midnside, nmaps, rcondlim)
          else
             call average_downgrade(map_in, map_mid, nside_in, midnside, nmaps)
          end if
       else
          map_mid(:,:) = map_in(:,:)
       end if

       if (smooth) then
          if (midnside /= nside_out) then
             call average_downgrade(map_mid, map_unsmoothed, midnside, nside_out, nmaps)
          else
             map_unsmoothed(:,:) = map_mid(:,:)
          end if

          IF (use_beam) THEN
             CALL windowfile_downgrade(file_beam, map_mid, map_out, midnside, nside_out, nmaps, .false., .true., .true.)
          ELSE IF (use_cosine_window) THEN
             CALL apodized_downgrade(map_mid, map_out, midnside, nside_out, nmaps, .false., .true., .true.)
          ELSE
             CALL smooth_downgrade(fwhm, map_mid, map_out, midnside, nside_out, nmaps, .false., .true., .true.)
          END IF

          if ( nmaps == 3 .and. .not. polsmooth ) then
             map_out(:,2:3) = map_unsmoothed(:,2:3)
          end if

          where (map_unsmoothed == 0) map_out = 0
       else
          map_out(:,:) = map_mid(:,:)
       end if
       

       DO imap = 1, nmaps
          evec_out((imap-1)*npix_out:imap*npix_out-1) = map_out(:,imap)
       END DO

       WRITE(evecUnit+1, rec=icol) evec_out(0:npix_out*nmaps-1)
    END DO

    CLOSE(unit=evecUnit)
    CLOSE(unit=evecUnit+1)

    deallocate(map_in, map_mid, map_out)

  END SUBROUTINE processEvecs



  SUBROUTINE parseArguments()

    IMPLICIT NONE
    
    LOGICAL :: there
    integer :: iargument

    IF (nArguments() < 2) THEN
       WRITE (*,'(/,a,/)') &
            '  Usage: downgrade_evecs <evec-file> <nside_out> ' // &
            '[-nobs <Nobs_matrices>] [-fwhm <FWHM>] [-beam <beamfile>] [-o <out>] ' // &
            '[ -nopolsmooth] [-rcond <rcond_lim>] [-midnside <midnside>]'
       STOP 
    END IF

    CALL getArgument(1, file_in)
    INQUIRE(file=TRIM(file_in), exist=there)
    IF (.NOT. there) STOP 'eigenvector file does not exist'
    WRITE (*,'(a)') ' Reading ' // trim(file_in)

    CALL getArgument(2, argument)
    READ (argument, *, iostat=ierr) nside_out
    if (ierr /= 0) then
       print *,'Failed to parse nside_out from ' // trim(argument)
       stop
    end if
    WRITE (*, '(a,i4)') ' Downgrading to nside_out == ', nside_out

    noiseWeight = .FALSE.
    smooth = .FALSE.
    polsmooth = .TRUE.
    midnside = -1
    use_cosine_window = .FALSE.
    file_out = '!out.fits'
    FWHM = 0
    rcondlim = 0
          
    iargument = 3
    DO
       IF (iargument > narguments()) EXIT
             
       CALL getArgument(iargument, argument)
       IF (INDEX(argument, '-nobs') /= 0) THEN
          iArgument = iArgument + 1           
          CALL getArgument(iargument, argument)
          noiseWeight = .TRUE.
          file_nobs = TRIM(argument)
          INQUIRE(file=TRIM(file_nobs), exist=there)
          IF (.NOT. there) STOP 'file_nobs does not exist'
          WRITE (*,'(a)') ' Doing noise weighted averaging according to ' // &
               TRIM(file_nobs)
       ELSE IF (INDEX(argument, '-beam') /= 0) THEN
          if (smooth) stop 'Cannot apply more then one smoothing window'
          smooth = .TRUE.
          use_beam = .TRUE.
          iArgument = iArgument + 1           
          CALL getArgument(iargument, argument)
          file_beam = TRIM(argument)
          INQUIRE(file=TRIM(file_beam), exist=there)
          IF (.NOT. there) STOP 'beamfile does not exist'
          WRITE (*,'(a)') ' Doing harmonic smoothing according to ' // TRIM(file_beam)
       ELSE IF (INDEX(argument, '-fwhm') /= 0) THEN
          iArgument = iArgument + 1           
          CALL getArgument(iargument, argument)
          READ (argument,*,iostat=ierr) FWHM
          IF (ierr /= 0) STOP 'Unable to parse FWHM'
          IF (FWHM > 1E-6) THEN
             if (smooth) stop 'Cannot apply more then one smoothing window'
             smooth = .TRUE.
             WRITE (*,'(a,g15.5,a)') &
                  ' Doing harmonic smoothing with FWHM == ', FWHM, ''' beam'
          ELSE IF (FWHM < -1E-6) THEN
             if (smooth) stop 'Cannot apply more then one smoothing window'
             smooth = .TRUE.
             use_cosine_window = .TRUE.
             WRITE (*,'(a,g15.5,a)') &
                  ' Doing harmonic smoothing with cosine window'
          ELSE
             WRITE (*,'(a)') ' FWHM=0 :  naive pixel averaging'
          END IF
       ELSE IF (INDEX(argument, '-rcond') /= 0) THEN
          iArgument = iArgument + 1           
          CALL getArgument(iargument, argument)
          READ (argument,*,iostat=ierr) rcondlim
          IF (ierr /= 0) STOP 'Unable to parse rcondlim'
          WRITE (*,'(a,es15.5)') ' Discarding pixels with rcond less than ', rcondlim
       ELSE IF (INDEX(argument, '-midnside') /= 0) THEN
          iArgument = iArgument + 1           
          CALL getArgument(iargument, argument)
          READ (argument,*,iostat=ierr) midnside
          IF (ierr /= 0) STOP 'Unable to parse midnside'
          WRITE (*,'(a,i8)') ' Using intermediate nside = ', midnside
       ELSE IF (INDEX(argument, '-o') /= 0) THEN
          iArgument = iArgument + 1           
          CALL getArgument(iargument, argument)
          file_out = TRIM(argument)
       ELSE IF (INDEX(argument, '-nopolsmooth') /= 0) THEN
          polsmooth = .FALSE.
          write (*,'(a)') ' ONLY smoothing the temperature component.'
       ELSE
          WRITE (*,*) 'Unrecognized argument: ' // TRIM(argument)
          stop
       END IF
       
       iargument = iargument + 1
    END DO
    
    if (noiseweight .and. (nside_in == nside_out)) &
         stop 'Cannot noise-weight without change in resolution'
    
    WRITE (*,'(a)') &
         ' Saving downgraded eigenvectors into '//TRIM(file_out)
 
  END SUBROUTINE parseArguments


END PROGRAM downgrade_evecs
