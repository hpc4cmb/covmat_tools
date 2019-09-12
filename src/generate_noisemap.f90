! This code reads in the eigenvectors and eigenvalues of a
! madam pixel-pixel covariance matrix.
! It then outputs a simulated noise map conforming to
! the eigenvectors and eigenvalues.
!
! February 10th 2009, RK

PROGRAM generate_noisemap

  USE healpix_types
  USE fitstools, ONLY   : output_map, getsize_fits
  USE head_fits, ONLY   : add_card, write_minimal_header
  USE pix_tools, ONLY   : convert_ring2nest
  USE rngmod
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, evec_file, eval_file
  CHARACTER(len=filenamelen) :: mapfile
  CHARACTER(len=80) :: header(1000)
  INTEGER(i4b) :: npix, nstokes, icol, nside
  INTEGER(i4b) :: seed
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: evec(:), evals(:), map(:, :), map2(:)
  REAL(dp) :: limit=1E20
  LOGICAL :: there, invert=.false.
  INTEGER(i4b) :: ierr, nhit, iargument
  TYPE(planck_rng) :: rnghandle
  REAL(dp) :: gauss

  IF (nArguments() < 5)  THEN
     WRITE (*,'(/,a,/)') &
          'Usage: generate_noisemap <covmat-evecs> <covmat-evals> ' // &
          '<nside> <nstokes> [-o <map>] [-lim <limit>] [-inv]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     evec_file = TRIM(argument)
     INQUIRE(file=TRIM(evec_file), exist=there)
     IF (.NOT. there) STOP 'evec_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(evec_file)

     CALL getArgument(2, argument)
     eval_file = TRIM(argument)
     INQUIRE(file=TRIM(eval_file), exist=there)
     IF (.NOT. there) STOP 'eval_file  does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(eval_file)

     CALL getArgument(3, argument)
     READ(argument, *) nside
     WRITE (*,*) ' nside = ', nside
     npix = 12*nside**2

     CALL getArgument(4, argument)
     READ(argument, *) nstokes
     WRITE (*,*) ' nstokes = ', nstokes

     iargument = 5
     limit = 0.0
     invert = .false.
     mapfile = '!noisemap.fits'
     DO
        IF (nArguments() < iargument) EXIT
        CALL getArgument(iargument, argument)

        IF (INDEX(argument, '-o') /= 0) THEN
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           mapfile = TRIM(ADJUSTL(argument))
        ELSE IF (INDEX(argument, '-lim') /= 0) THEN
           iArgument = iArgument + 1
           IF (nArguments() < iargument) STOP 'Must supply argument with -lim'
           CALL getArgument(iArgument, argument)
           READ(argument, *, iostat=ierr) limit
           IF (ierr /= 0) stop 'Cannot parse limit'
           WRITE (*,*) 'Using threshold eval > evalmax * ',limit
        ELSE IF (INDEX(argument, '-inv') /= 0) THEN
           invert = .TRUE.
           WRITE (*,*) 'Inverting eigenvalues'
        ELSE
           WRITE (*,*) 'Unrecognized command line option: ',TRIM(argument)
           STOP
        END IF

        iargument = iargument + 1
     END DO

     WRITE (*,*) ' Generating ' // TRIM(mapfile)
  END IF

  ALLOCATE(map(0:npix-1,nstokes), map2(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for maps'

  CALL tic
  ALLOCATE(evec(0:npix*nstokes-1), evals(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for evecs'

  OPEN(unit=covmat_unit, file=TRIM(eval_file), status='old', &
       form='formatted')
  READ(covmat_unit,*) evals
  CLOSE(covmat_unit)

  WHERE(evals < maxval(evals)*limit) evals = 0.0_dp
  IF (invert) THEN
     ! assume that evals correspond to N^-1
     WHERE(evals /= 0.0_dp) evals= 1.0/evals
  END IF

  OPEN(unit=covmat_unit, file=TRIM(evec_file), status='old', &
       form='unformatted', access='direct', recl=npix*nstokes*8)
  nhit = 0

  CALL SYSTEM_CLOCK(count=seed)
  CALL rand_init(rnghandle, seed)
  map2 = 0.0
  DO icol = 0, npix*nstokes - 1
     IF (evals(npix*nstokes-1) > 0.0_dp) THEN
        READ(covmat_unit, rec=icol+1) evec
        gauss = rand_gauss(rnghandle)
        nhit  = nhit + 1
        map2  = map2 + gauss*evals(icol)**0.5*evec
     END IF
  END DO
  
  IF (nhit < npix*nstokes) WRITE (*,*) 'Excluded ', npix*nstokes-nhit, &
       ' eigenvectors.'

  CALL toc('Generate map')

  map = RESHAPE(map2, (/npix, nstokes/))

  CALL write_minimal_header(header, 'MAP', nside=nside, ordering='NESTED', &
       creator='generate_noisemap', randseed=seed, polar=(nstokes>1))

  CALL output_map(map, header, mapfile)

  DEALLOCATE(map, map2, evec, evals)

END PROGRAM generate_noisemap
