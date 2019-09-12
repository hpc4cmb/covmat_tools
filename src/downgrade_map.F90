PROGRAM downgrade_map

  USE fitstools
  USE head_fits
  USE pix_tools
  USE healpix_types
  USE extension

  USE covmat_util, ONLY : tic, toc
  USE downgrade

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: filename_in, filename_out, file_nobs, file_beam
  INTEGER(i4b) :: nside_in, nMaps, nside_out, midnside=-1
  INTEGER(i8b) :: npix_in, npix_out
  REAL(dp) :: FWHM, apodize
  REAL(dp), POINTER :: map_in(:, :), map_out(:, :), map_unsmoothed(:, :)
  REAL(sp), POINTER :: map_out_sp(:, :)
  CHARACTER(len=80) :: header(1000), kwd, val
  INTEGER(i4b) :: ordering, ierr, iMap, istart
  INTEGER(i8b) :: npixtot, npix2one
  LOGICAL :: there, noiseWeight=.FALSE., output_double_precision=.FALSE., use_beam=.FALSE.
  LOGICAL :: smooth=.FALSE., use_cosine_window=.FALSE., do_hits=.FALSE.
  LOGICAL :: polsmooth=.TRUE., invert=.FALSE., diag_only=.false.
  REAL(dp) :: rcondlim=0

  CALL tic(263)

  CALL parseArguments()

  IF (do_hits .AND. (noiseWeight .OR. smooth)) &
       STOP 'Cannot noise-weight or smooth hits'

  CALL tic
  CALL readInputMap(map_in)
  CALL toc('Read map')

  if (apodize > 0) then
     CALL tic
     CALL apodize_map(map_in, npix_in, nside_in, nmaps, apodize)
     CALL toc('Apodize map')
  end if

  IF (smooth .and. (nmaps /= 1 .and. nmaps /= 3)) stop 'Do not know how to smooth nmaps /= 1,3'

  IF (noiseweight .AND. (nside_in < nside_out)) &
       STOP 'Cannot noise-weight to higher resolution'

  !call subtract_monopole(map_in)

  ! if we are both noise-weighting and smoothing, choose an intermediate
  ! resolution to reduce aliasing effects
  if (smooth) then
     if (midnside > 0) then
        if (midnside > nside_in) stop 'midnside must be less than input nside'
        if (midnside < nside_out) stop 'midnside must be greater than or equal the output nside'
     else
        IF (nside_out < nside_in) THEN
           midnside = nside_out * 2
        ELSE
           midnside = nside_out
        END IF
     end if
  else
     if (midnside > 0) print *,'Warning, midnside ignored without smoothing'
     midnside = nside_out
  end if

  npix_out = nside2npix(midnside)
  npix2one = npix_in / npix_out

  CALL tic
  IF (midnside < nside_in) THEN
     ALLOCATE(map_out(0:npix_out-1, nMaps), stat=ierr)
     IF (ierr /= 0) STOP 'No room for output map !?'
     WRITE (*,'(/,a,i5,a,/)') ' Downgrading ', npix2one, ' pixels into one.'

     IF (noiseweight) THEN
        WRITE (*,'("Noise-weighting to nside == ",i0)') midnside
        CALL noise_weight_downgrade(file_nobs, map_in, map_out, &
             nside_in, midnside, nmaps, rcondlim, invert, diag_only)
     ELSE
        WRITE (*,'("Averaging to nside == ",i0)') midnside
        CALL average_downgrade(map_in, map_out, &
             nside_in, midnside, nmaps)
        IF (do_hits) map_out = INT(map_out * npix2one, i4b) ! from average to sum
     END IF
     CALL toc('Downgrade')

     IF (smooth) THEN
        CALL tic
        DEALLOCATE(map_in)
        ALLOCATE(map_in(0:npix_out-1, nMaps), stat=ierr)
        IF (ierr /= 0) STOP 'No room for output map !?'
        map_in = real(map_out, sp)
        nside_in = midnside
        npix_in = nside2npix(midnside)
        DEALLOCATE(map_out)
     END IF
  END IF

  IF (smooth) THEN

     if (apodize > 0) then
        CALL tic
        CALL apodize_map(map_in, npix_in, nside_in, nmaps, apodize)
        CALL toc('Apodize map')
     end if

     npix_out = nside2npix(nside_out)
     npix2one = npix_in/npix_out

     ! Produce downgraded version without smoothing

     ALLOCATE(map_unsmoothed(0:npix_out-1, nMaps), stat=ierr)
     IF (ierr /= 0) STOP 'No room for unsmoothed map !?'
     IF (nside_out == midnside) THEN
        map_unsmoothed(:, :) = map_in(:, :)
     ELSE
        ! Use average downgrading for the last step. This is a better match to
        ! what happens to the signal in the smoothing case
        CALL average_downgrade(map_in, map_unsmoothed, nside_in, nside_out, nmaps)
        CALL toc('Downgrade again')
     END IF

     ! Produce downgraded version with smoothing

     WRITE (*,'("Smoothing to nside== ",i0)') nside_out
     ALLOCATE(map_out(0:npix_out-1, nMaps), stat=ierr)
     IF (ierr /= 0) STOP 'No room for output map !?'
     WRITE (*,'(/,a,i5,a,/)') ' Downgrading ', npix2one, ' pixels into one.'

     WRITE (*,'(" Smoothing the map..")')
     if (use_beam) then
        CALL windowfile_downgrade(file_beam, map_in, map_out, nside_in, &
             nside_out, nmaps)
        CALL toc('Beam smooth')
     else IF (use_cosine_window) THEN
        CALL apodized_downgrade(map_in, map_out, nside_in, nside_out, nmaps)
        CALL toc('Apodized smooth')
     ELSE
        CALL smooth_downgrade(FWHM, map_in, map_out, nside_in, nside_out, nmaps)
        CALL toc('Gaussian smooth')
     END IF

     IF (nmaps == 3 .AND. .NOT. polsmooth) THEN
        map_out(:, 2:3) = map_unsmoothed(:, 2:3)
     END IF

     WHERE (map_unsmoothed == 0) map_out = 0

  END IF
  CALL toc('downgrade map')

  PRINT *,'nside_out : ' ,nside_out
  PRINT *,'npix_out : ', npix_out

  CALL tic
  header = ''
  call write_minimal_header(header, 'map', nside=nside_out, order=2, &
       fwhm_degree=fwhm/60.0, polar=(nmaps==3), creator='downgrade_map')
  CALL add_card(header, 'ORIG_NSIDE', nside_in, 'Resolution before downgrading')
  IF (noiseWeight) then
     ! Have to cut the absolute path from the name of the nobs file:
     ! the HEASARC convention for CONTINUE lines is not part of the standard and produces problems
     istart = index(file_nobs, '/', back=.true.)
     if (istart == 0) istart = 1
     CALL add_card(header, 'FILENOBS', TRIM(file_nobs(istart:)), '3x3 Nobs file')
  end IF
  if (do_hits) then
     call add_card(header, 'TTYPE1', 'HITS')
  else
     if (nmaps == 1 .or. nmaps == 3) then
        call add_card(header, 'TTYPE1', 'I_STOKES')
        if (nmaps == 3) then
           call add_card(header, 'TTYPE2', 'Q_STOKES')
           call add_card(header, 'TTYPE3', 'U_STOKES')
        end if
     else
        do imap = 1, nmaps
           write(kwd,'("TTYPE",i0)') imap
           write(val,'("COMP",i0)') imap
           call add_card(header, kwd, val)
        end do
     end if
  end if

  IF (output_double_precision) THEN
     WHERE(map_out == 0) map_out = HPX_DBADVAL
     CALL output_map(map_out, header, filename_out)
  ELSE
     ALLOCATE(map_out_sp(0:npix_out-1, nMaps), stat=ierr)
     map_out_sp = real(map_out, sp)
     WHERE(map_out_sp == 0) map_out_sp = HPX_SBADVAL
     CALL output_map(map_out_sp, header, filename_out)
     DEALLOCATE(map_out_sp)
  END IF
  CALL toc('write map')

  CALL toc('total execution', 263)

  DEALLOCATE(map_out)

CONTAINS



!!$  SUBROUTINE subtract_monopole(map)
!!$    IMPLICIT NONE
!!$
!!$    REAL(dp), POINTER :: map (:,:)
!!$    REAL(dp)          :: monopole
!!$
!!$    monopole = SUM(map(:,1))/(UBOUND(map,1)+1)
!!$    map(:,1) = map(:,1) - monopole
!!$
!!$    WRITE(*,*) 'Subtracted I-monopole: ', monopole
!!$  END SUBROUTINE subtract_monopole



  SUBROUTINE readInputMap(map_in)

    IMPLICIT NONE

    REAL(dp), POINTER :: map_in(:, :)

    integer :: nmaps_tot

    header = ''

    npixtot = getsize_fits(filename_in, nMaps=nMaps_tot, ordering=ordering, &
         nside=nside_in)
    npix_in = nside2npix(nside_in)
    if (nmaps < 1) nmaps = nmaps_tot

    WRITE (*,'(a,i9)')  ' nside_in == ', nside_in
    WRITE (*,'(a,i9)')  ' npix_in  == ', npix_in
    WRITE (*,'(a,i9)')  ' nMaps    == ', nMaps
    WRITE (*,'(a,i9)')  ' ordering == ', ordering

    IF (nside_in < nside_out) STOP 'Cannot downgrade to higher resolution'

    npix_in  = nside2npix(nside_in)

    ALLOCATE(map_in(0:npix_in-1, nMaps), stat=ierr)
    IF (ierr /= 0) STOP 'No room to store the maps!?'

    if (nmaps_tot == 2 .and. nmaps == 3) then
       ! Special case, were are downgrading a polarization correction map
       map_in(:,1) = 0
       CALL input_map(filename_in, map_in(:, 2:3), npix_in, 2, fmissval=0d0, &
            header=header)
    else
       CALL input_map(filename_in, map_in, npix_in, nMaps, fmissval=0d0, &
            header=header)
    end if

    SELECT CASE (ordering)
    CASE (0)
       WRITE (*,*) 'Unknown ordering, assuming NESTED'
    CASE (1)
       WRITE (*,*) 'Converting input map to NESTED'
       CALL convert_ring2nest(nside_in, map_in)
       CALL add_card(header, 'ORDERING', 'NESTED', 'Pixel ordering scheme')
    CASE (2)
       WRITE (*,*) 'Input map is in NESTED scheme'
    END SELECT

  END SUBROUTINE readInputMap


  subroutine apodize_map(map, npix, nside, nmaps, apodlim)
    integer(i8b) :: npix
    integer(i4b) :: nside, nmaps
    REAL(dp) :: map(0:npix-1, nmaps)
    real(dp) :: apodlim

    real(dp) :: mean, var, sigma
    integer :: imap, ierr, ngood, nlist, n_apod
    integer(i8b) :: pix, list(8)
    logical, allocatable :: mask(:)

    stop 'Apodization details still open'

    allocate(mask(0:npix-1), stat=ierr)
    if (ierr /= 0) stop 'No room for apodizing'

    mask = map(:, 1) /= 0
    ngood = count(mask)

    do imap = 1, nmaps
       n_apod = 0
       mean = sum(map(:, imap), mask=mask) / ngood
       var = sum((map(:, imap)-mean)**2, mask=mask) / ngood
       sigma = sqrt(var)
       do pix = 0, npix-1
          if (.not. mask(pix)) cycle
          if (abs(map(pix,imap)-mean) > apodlim*sigma) then
             CALL neighbours_nest(nside, pix, list, nlist)
             map(pix, imap) = sum(map(list(1:nlist), imap)) / nlist
             n_apod = n_apod + 1
          end if
       end do
       print *,' Apodized ', n_apod, ' pixels in map # ', imap
    end do

    deallocate(mask)

  end subroutine apodize_map



  SUBROUTINE parseArguments()

    IMPLICIT NONE

    INTEGER :: iargument, ierr
    CHARACTER(filenamelen) :: argument

    IF (nArguments() < 3) THEN
       WRITE (*,'(/,a,/)') '  Usage: downgrade_map <map> <nside_out> ' // &
            '[-nobs <Nobs_matrices>] ' // &
            '[-wcov <white_noise_covariace_matrices] ' // &
            '[-diagnobs <Nobs_matrices>] ' // &
            '[-diagwcov <white_noise_covariace_matrices] ' // &
            '[-fwhm <FWHM>] [-beam <beamfile>] [-o <out>] [-d] [--hits] ' // &
            '[-nopolsmooth] [-rcondlim <rcondlim>] [-midnside <midnside>] ' // &
            '[-nmaps <nmaps>] [-apodize <sigma>]'
       STOP
    END IF

    CALL getArgument(1, filename_in)
    INQUIRE(file=TRIM(filename_in), exist=there)
    IF (.NOT. there) STOP 'map does not exist'

    CALL getArgument(2, argument)
    READ (argument, *, iostat=ierr) nside_out
    IF (ierr /= 0) THEN
       PRINT *,'Failed to parse nside_out from ' // TRIM(argument)
       STOP
    END IF
    WRITE (*, '(a,i4)') ' Downgrading to nside_out == ', nside_out

    noiseWeight = .FALSE.
    smooth = .FALSE.
    use_cosine_window = .FALSE.
    use_beam = .false.
    invert = .FALSE.
    diag_only = .false.
    filename_out = '!out.fits'
    rcondlim = 0
    nmaps = 0
    apodize = 0

    iargument = 3
    FWHM = 0
    DO
       IF (iargument > narguments()) EXIT

       CALL getArgument(iargument, argument)
       IF (INDEX(argument, 'nobs') /= 0) THEN
          diag_only = INDEX(argument, 'diag') /= 0
          iArgument = iArgument + 1
          CALL getArgument(iargument, argument)
          noiseWeight = .TRUE.
          file_nobs = TRIM(argument)
          INQUIRE(file=TRIM(file_nobs), exist=there)
          IF (.NOT. there) STOP 'file_nobs does not exist'
          if (diag_only) then
             WRITE (*,'(a)') ' Doing noise weighted averaging according to ' &
                  // 'diagonal of ' // TRIM(file_nobs)
          else
             WRITE (*,'(a)') ' Doing noise weighted averaging according to ' &
                  // TRIM(file_nobs)
          end if
       ELSE IF (INDEX(argument, 'beam') /= 0) THEN
          if (smooth) stop 'Cannot apply more then one smoothing window'
          smooth = .TRUE.
          use_beam = .TRUE.
          iArgument = iArgument + 1
          CALL getArgument(iargument, argument)
          file_beam = TRIM(argument)
          INQUIRE(file=TRIM(file_beam), exist=there)
          IF (.NOT. there) STOP 'beamfile does not exist'
          WRITE (*,'(a)') ' Doing harmonic smoothing according to ' // TRIM(file_beam)
       ELSE IF (INDEX(argument, 'wcov') /= 0) THEN
          diag_only = INDEX(argument, 'diag') /= 0
          invert = .TRUE.
          iArgument = iArgument + 1
          CALL getArgument(iargument, argument)
          noiseWeight = .TRUE.
          file_nobs = TRIM(argument)
          INQUIRE(file=TRIM(file_nobs), exist=there)
          IF (.NOT. there) STOP 'file_wcov does not exist'
          if (diag_only) then
             WRITE (*,'(a)') ' Doing noise weighted averaging according to ' &
                  // 'the diagonal of inverse of ' // TRIM(file_nobs)
          else
             WRITE (*,'(a)') ' Doing noise weighted averaging according to ' &
                  // 'the inverse of ' // TRIM(file_nobs)
          end if
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
             smooth          = .TRUE.
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
          WRITE (*,'(a,es15.5)') ' Discarding pixels with rcond less than ',rcondlim
       ELSE IF (INDEX(argument, '-midnside') /= 0) THEN
          iArgument = iArgument + 1
          CALL getArgument(iargument, argument)
          READ (argument,*,iostat=ierr) midnside
          IF (ierr /= 0) STOP 'Unable to parse midnside'
          WRITE (*,'(a,i8)') ' Using intermediate nside = ', midnside
       ELSE IF (INDEX(argument, '-hits') /= 0) THEN
          do_hits = .TRUE.
          WRITE (*,*) 'Map is a hit map'
       ELSE IF (INDEX(argument, '-o') /= 0) THEN
          iArgument = iArgument + 1
          CALL getArgument(iargument, argument)
          filename_out = TRIM(argument)
       ELSE IF (INDEX(argument, '-d') /= 0) THEN
          output_double_precision = .TRUE.
          WRITE (*,'(a)') ' Output map will be in DOUBLE precision.'
       ELSE IF (INDEX(argument, '-nopolsmooth') /= 0) THEN
          polsmooth = .FALSE.
          WRITE (*,'(a)') ' ONLY smoothing the temperature component.'
       ELSE IF (INDEX(argument, '-nmap') /= 0) THEN
          iArgument = iArgument + 1
          CALL getArgument(iargument, argument)
          READ (argument,*,iostat=ierr) nmaps
          IF (ierr /= 0) STOP 'Unable to parse nmaps'
          WRITE (*,'(a,i0,a)') ' Will read first ',nmaps,' maps'
       ELSE IF (INDEX(argument, '-apodize') /= 0) THEN
          iArgument = iArgument + 1
          CALL getArgument(iargument, argument)
          READ (argument,*,iostat=ierr) apodize
          IF (ierr /= 0) STOP 'Unable to parse apodize'
          WRITE (*,'(a,es15.5,a)') ' Will apodize values beyond ',apodize,' sigma'
       ELSE
          WRITE (*,*) 'Unrecognized argument: ' // TRIM(argument)
          STOP
       END IF

       iargument = iargument + 1
    END DO

    WRITE (*,'(a)') ' Saving downgraded map into ' // TRIM(filename_out)

  END SUBROUTINE parseArguments


END PROGRAM downgrade_map
