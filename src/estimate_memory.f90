! This small tool estimates the memory needed to compute
! the inverse noise covariance from Madam

PROGRAM estimate_memory

  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc
  USE fourier

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: key, value, covmat_file, outfile
  INTEGER(I8B) :: ndays, nstokes=3, nsamples_per_baseline, ndets=2
  INTEGER(I8B) :: nsamples, nbaselines, nside, nprocs
  INTEGER(I8B) :: npix_per_ring, npix, ndiagonal, mem_per_proc
  REAL(dp) :: fsample, fknee, fmin=0.15741E-5, alpha=1.7
  REAL(dp) :: memsize, rtemp, pixelsize
  REAL(dp) :: zero_limit=1E-8
  INTEGER(I4B) :: ierr, iArgument

  IF (nArguments() < 2)  THEN
     WRITE (*,*) 'Usage: estimate_memory ' // &
          '--fsample <fsample> --fknee <fknee> --days <days> ' // &
          '--bllength <bllength> --nside <nside> [--nstokes <nstokes>]'
     stop
  ELSE
     iArgument = 1
     DO
        CALL getArgument(iArgument,   key)
        CALL getArgument(iArgument+1, value)

        if (index(key, '--fsample') /= 0) then
           read (value,*) fsample
           write (*,'(a,g10.3,a)') '  fsample == ', fsample, ' Hz'
        else if (index(key, '--fknee') /= 0) then
           read (value,*) fknee
           write (*,'(a,ES10.3,a)') '    fknee == ', fknee, ' Hz'
        else if (index(key, '--alpha') /= 0) then
           read (value,*) alpha
           write (*,'(a,ES10.3)')   '    alpha == ', alpha
        else if (index(key, '--fmin') /= 0) then
           read (value,*) fmin
           write (*,'(a,ES10.3,a)') '     fmin == ', fmin, ' Hz'
        else if (index(key, '--nside') /= 0) then
           read (value,*) nside
           write (*,'(a,i10)')    '    nside == ', nside
           npix = 12*nside**2
        else if (index(key, '--days') /= 0) then
           read (value,*) ndays
           write (*,'(a,i10)')    '     days == ', ndays
        else if (index(key, '--bllength') /= 0) then
           read (value,*) nsamples_per_baseline
           write (*,'(a,i10,a)')    ' bllength == ', &
                nsamples_per_baseline, ' samples'
        else
           write (*,*) 'unrecognized key: ' // trim(key)
           stop
        end if
        iArgument = iArgument + 2
        if (iArgument >  nArguments() - 1) exit
     END DO
  END IF

  memsize = 0
  ! perform the approximations
  nsamples = ndays*24*3600*fsample
  nbaselines = nsamples / nsamples_per_baseline

  rtemp = nsamples*3.0*8/2**20
  write (*,'(/,a,i8,a)') ' Memory for pointing : ', int(rtemp, i8b), 'MB'
  memsize = memsize + rtemp

  pixelsize = sqrt(FOURPI/npix)
  npix_per_ring = TWOPI / pixelsize
  write (*,'(a,ES10.3,a)') ' Pixel size is ', RAD2DEG*pixelsize, ' deg'
  write (*,'(a,I10,a)')    ' translates to ', npix_per_ring, &
       ' pixels per scanning ring'
  write (*,*) ' (== hit_pixels_per_detector, when nprocs > npix_per_ring)'

  ndiagonal = get_diagonal_width(fsample, fknee, fmin, alpha, &
       nsamples_per_baseline, zero_limit)

  write (*,*) 'Diagonal width is ', ndiagonal

  nprocs = 1
  mem_per_proc = 0
  write (*,'(/,a,i2,a/)') &
       ' Memory requirements per processor (', ndets, ' detectors): '
  do
     mem_per_proc = npix_per_ring * nstokes * (nbaselines/nprocs + 2 * ndiagonal) * ndets ! PTFA
     mem_per_proc = mem_per_proc + npix_per_ring * nstokes * nbaselines/nprocs * ndets   ! PTF
     mem_per_proc = mem_per_proc + nsamples * nstokes / nprocs ! pointing
     mem_per_proc = mem_per_proc * 8 / 2**20
     write (*,'(i12,a,i12,a)') nprocs, ' : ', mem_per_proc, 'MB'
     if (mem_per_proc < 100) exit
     nprocs = nprocs * 2
  end do



contains



  function get_diagonal_width(fsample, fknee, fmin, alpha, &
       nsamples_per_baseline, zero_limit)
    ! This routine uses the instrument noise model to study the middlematrix
    ! diagonal width
    real(dp), intent(in)     :: fsample, fknee, fmin, alpha, zero_limit
    integer(i8b), intent(in) :: nsamples_per_baseline
    integer(i8b)             :: get_diagonal_width
    real(dp)                 :: sigma=1.0, fa, x
    integer(i8b), parameter  :: length = 2**20
    complex(dpc), pointer    :: psd(:) ! noise power spectral density
    integer(i8b)             :: i, k
    integer(i4b)             :: ierr, slength=length
    real(dp), pointer        :: aspec(:), f(:), g(:), spectrum(:)
    real(dp), pointer        :: ftable(:),spectrum_table(:,:)
    real(dp), pointer        :: cov_a_inv(:)

    allocate(psd(length/2+1), aspec(length/2+1), spectrum(length), &
         f(length), g(length), cov_a_inv(length), stat=ierr)
    if (ierr /= 0) stop 'no room for psd'
    
    ! initialize psd using the detector noise model
    ! this part is from Madam
    fa    = fsample / nsamples_per_baseline
    aspec = 0.0
    do k = 0,10
       do i = 1, length
          f(i) = k*fa + (i-1)*fa/length
       end do
       
       do i = 1, length
          spectrum(i) = sigma**2/fsample*&
               fknee**alpha/(f(i)**alpha+fmin**alpha)
       end do
       
       do i = 1, length
          x = pi*f(i)/fa
          if (x.lt.1.e-30) then
             g(i) = 1.d0
          else
             g(i) = sin(x)**2/x**2
          end if
       end do

       if (k == 0) then
          aspec(1) = spectrum(1)
       else
          aspec(1) = aspec(1) + 2*spectrum(1)*g(1)
       end if

       do i = 2, length/2+1
          aspec(i) = aspec(i) + spectrum(i)*g(i) + &
               spectrum(length-i+2)*g(length-i+2)
       end do
    end do

    psd  = 1.0/(fa*aspec)

    CALL init_fourier(slength)
    CALL dfftinv(cov_a_inv, psd)

    ! Invert the middlematrix by fourier transform
    ! store results in cov_a_inv
    cov_a_inv(1) = cov_a_inv(1) + nsamples_per_baseline
    call dfft(psd, cov_a_inv)
    call dfftinv(cov_a_inv, 1/psd)

    !write (*,*) cov_a_inv(1:5)
    !write (*,*) cov_a_inv(length-5:length)

    call close_fourier()

    do get_diagonal_width = 2, length
       if (cov_a_inv(get_diagonal_width)/cov_a_inv(1) < zero_limit) return
    end do

  end function get_diagonal_width

END PROGRAM estimate_memory
