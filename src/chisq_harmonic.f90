! This code reads madam alm-alm inverse covariance matrix and a
! given residual noise map. It then computes the chi squared
! test for the given map using the inverse covariance matrix
!
! Matrix is assumed not to include ell<2 and m<0 modes
!
! 2008-08-27 RK

PROGRAM chisq_test

  USE healpix_types
  USE fitstools, ONLY   : input_map, getsize_fits
  USE pix_tools, ONLY   : convert_nest2ring
  USE alm_tools, ONLY   : map2alm
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, mapfile, form
  INTEGER(i4b) :: npix, nstokes, icol, nside, ordering
  INTEGER(i4b) :: nalm, nalm2, lmin, lmax, mmax, m, ell, count
  INTEGER(i8b) :: dnpix
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: map(:, :), w8ring(:, :)
  COMPLEX(dpc), POINTER :: covmat(:), chisq(:, :), almIn(:, :, :)
  COMPLEX(dpc), POINTER :: alm(:, :), alm2(:, :, :), alm3(:), alm4(:)
  REAL(dp) :: dof, zbounds(2), chisq_tot
  LOGICAL :: there
  INTEGER(i4b) :: ierr, isig1, isig2, isig, rec_len

  IF (nArguments() /= 4)  THEN
     STOP 'Usage: chisq_harmonic <inv_covmat> <lmin> <lmax> <map>'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) THEN
        WRITE (*,*) TRIM(covmat_file) // ' does not exist!'
        STOP 
     END IF
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) lmin
     WRITE (*,'(a,i4)') '  lMin    == ', lmin
     mmax = lmax

     CALL getArgument(3, argument)
     READ (argument, *) lmax
     WRITE (*,'(a,i4)') '  lMax    == ', lmax
     mmax = lmax

     CALL getArgument(4, argument)
     mapfile = TRIM(argument)
     INQUIRE(file=mapfile, exist=there)
     IF (.NOT. there) THEN
        WRITE (*,*) TRIM(mapfile) // ' does not exist!'
        STOP 
     END IF
     WRITE (*,*) ' Reading ' // TRIM(mapfile)
  END IF

  ! Read in the map
  CALL tic
  dnpix = getsize_fits(mapfile, nmaps=nstokes, ordering=ordering, &
       nside=nside)
  npix = int(dnpix, i4b)
  ALLOCATE(map(0:npix-1,nstokes), stat=ierr)
  IF (ierr /= 0) STOP 'No room for map'
  CALL input_map(mapfile, map, npix, nstokes)
  IF (ordering == 2) THEN
     ! map2alm assumes map to be in RING scheme
     WRITE (*,*) ' Note: converting map into RING scheme'
     CALL convert_nest2ring(nside, map)
  END IF
  CALL toc('read map')


  ! find the harmonic expansion
  nalm  = (lmax+1)**2 - lmin**2 ! all alms
  nalm2 = nalm - (nalm - (lmax+1-lmin))/2 ! exclude degenerate modes
  WRITE (*,'(a,i6)') 'npix    == ', npix
  WRITE (*,'(a,i6)') 'nalm(lmin,lmax)     == ', nalm
  WRITE (*,'(a,i6)') 'nalm(nondegenerate) == ', nalm2
  ALLOCATE(almIn(nstokes, 0:lmax, 0:mmax), alm(nalm2, nstokes), &
       alm2(nalm2, nstokes, nstokes), w8ring(2*nside, nstokes), stat=ierr)
  IF (ierr /= 0) STOP 'no room for alm'
  zbounds = 0.0
  w8ring  = 1.0
  CALL map2alm(nside, lmax, mmax, map, almIn, zbounds, w8ring)
  DO isig = 1, nstokes
     count = 1
     DO ell = lmin, lmax
        DO m = 0, ell
           alm(count, isig) = almIn(isig, ell, m)
           count            = count + 1
        END DO
     END DO
  END DO

  ! As a test generate alms with exactly desired covariance:
  !call generate_alms(alm,nalm2,nstokes,'eigenvalues.dat','eigenvectors.dat')


  ! Perform (m^T N^-1 m) by storing only a single column of
  ! noise covariance at a time. To report a breakdown of the
  ! chisq, this is done in blocks
  CALL tic
  ALLOCATE(covmat(nalm2*nstokes), chisq(nstokes, nstokes), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  DO icol = 1, nalm2*nstokes
     READ(covmat_unit, rec=icol) covmat
     isig1 = (icol-1)/nalm2 + 1
     DO isig2 = 1, nstokes
        alm2(MODULO(icol-1,nalm2)+1, isig2, isig1) = &
             DOT_PRODUCT(alm(:,isig2), covmat((isig2-1)*nalm2+1:isig2*nalm2))
     END DO
  END DO
  DO isig1 = 1, nstokes
     DO isig2 = 1, nstokes
        chisq(isig2, isig1) = &
             DOT_PRODUCT(CONJG(alm2(:, isig2, isig1)), alm(:, isig1))
     END DO
  END DO
  CLOSE(covmat_unit)
  CALL toc('compute chi-squared')

  form = '(a,10es11.2)'
  WRITE (*,*)
  WRITE (*,form) '             TT TE TB      ', chisq(1,1:3)
  WRITE (*,form) ' chisq   ==  ET EE EB  ==  ', chisq(2,1:3)
  WRITE (*,form) '             BT BE BB      ', chisq(3,1:3)
  WRITE (*,*)

  chisq_tot = REAL(SUM(chisq), dp)
  WRITE (*,form) ' Total chi squared ==  ', chisq_tot

  WRITE (*,*)
  dof = nalm2*nstokes
  WRITE (*,form) ' chi squared / dof ==  ', chisq_tot/dof
  WRITE (*,*)
  
  WRITE (*,form) ' Expectance == ', REAL(dof,dp)
  WRITE (*,form) '      stdev == ', (2*dof)**0.5
  WRITE (*,'(/,a,G15.5,a)') &
       ' Deviation from chi squared expected value : ', &
       REAL((chisq_tot-dof),dp)/(2*dof)**0.5, ' sigma'
  WRITE (*,*)

  ALLOCATE(alm3(nstokes*nalm2), alm4(nstokes*nalm2))
  DO isig = 1, nstokes
     alm3((isig-1)*nalm2+1:isig*nalm2) = alm(:,isig)
  END DO
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  DO icol = 1, nalm2*nstokes
     READ(covmat_unit, rec=icol) covmat
     alm4(icol) = DOT_PRODUCT(alm3, covmat)
  END DO
  CLOSE(covmat_unit)
  WRITE (*,*) DOT_PRODUCT(CONJG(alm4), alm3)
  

!!$CONTAINS
!!$
!!$
!!$  SUBROUTINE generate_alms(alms, nalm, nstokes, evalue_file, evec_file)
!!$    ! Use the eigenvectors of the covariance matrix convert independet
!!$    ! random variables into having the desired covariance properties.
!!$    ! In effect the generated alm-realization is a linear combination
!!$    ! of the eigenvectors, where the weights are sqrt(lambda)*rand
!!$    USE rngmod, ONLY : rand_init, rand_gauss, planck_rng
!!$
!!$    IMPLICIT NONE
!!$
!!$    COMPLEX(dpc), POINTER   :: alms(:, :)
!!$    INTEGER(i4b) :: nalm, nstokes
!!$    CHARACTER(len=*) :: evalue_file, evec_file
!!$    INTEGER(i4b), PARAMETER :: evalue_unit=67, evec_unit=68
!!$    COMPLEX(dpc), POINTER :: evec(:)
!!$    TYPE(planck_rng), SAVE :: rng_handle
!!$    LOGICAL, SAVE :: initialized = .FALSE.
!!$    REAL(dp) :: gauss1, gauss2, evalue
!!$    INTEGER(i4b) :: count, rate, icol, ierr, idummy
!!$    INTEGER(i4b) :: rec_len
!!$
!!$    IF (.NOT. initialized) THEN
!!$       CALL SYSTEM_CLOCK(count, rate)
!!$       CALL rand_init(rng_handle, count)
!!$       initialized = .TRUE.
!!$    END IF
!!$
!!$    ALLOCATE(evec(nalm), stat=ierr)
!!$    IF (ierr /= 0) STOP 'No room for evec'
!!$
!!$    OPEN(unit=evalue_unit, file=evalue_file, status='old', form='formatted')
!!$    
!!$    inquire(iolength=rec_len) evec
!!$    OPEN(unit=evec_unit, file=evec_file, status='old', form='unformatted', &
!!$         access='direct', recl=rec_len)
!!$    alms = 0
!!$    count = 1
!!$    DO icol = 1, nalm*nstokes
!!$       READ(evalue_unit, *) idummy, evalue
!!$       gauss1 = rand_gauss(rng_handle)/2**0.5
!!$       gauss2 = rand_gauss(rng_handle)/2**0.5
!!$       DO isig = 1, nstokes
!!$          READ(evec_unit, rec=count) evec
!!$          !write (*,'(f7.3)', advance='no'), abs(dot_product(evec, evec))
!!$          alms(:,isig) = alms(:,isig) &
!!$               + evalue**0.5*CMPLX(gauss1, gauss2, dpc)*evec
!!$          count = count + 1
!!$        END DO
!!$        !write (*,*)
!!$    END DO
!!$  END SUBROUTINE generate_alms
       
END PROGRAM chisq_test
