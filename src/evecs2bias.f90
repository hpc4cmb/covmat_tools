PROGRAM evecs2bias

!!$ 14 November, 2007 / Torsti Poutanen
!!$ Modified for cmbtools and large files 4 July 2008 / Reijo Keskitalo

  USE healpix_types
  USE pix_tools
  USE fitstools
  use head_fits
  USE alm_tools
  USE extension, ONLY : getArgument, nArguments

  IMPLICIT NONE

  INTEGER(i4b) :: Npix, NNpix, nside, lmax
  INTEGER(i4b) :: i, l, evec_unit=543, iargument, imap
  real(sp), allocatable :: mask(:, :)
  REAL(dp), ALLOCATABLE :: dmapout(:, :), nl(:,:)
  REAL(dp), ALLOCATABLE :: eigenvalue(:), nltest(:,:), w8ring_TQU(:, :)
  REAL(dp), ALLOCATABLE :: eigenvector(:)
  COMPLEX(dpc), ALLOCATABLE :: alm_TGC_test(:, :, :)
  REAL(dp) :: zbounds(1:2), limit=0
  REAL(sp) :: t0, t3, t4
  INTEGER(i4b) :: nmaps, nlheader, ncl, nPlm, ierr, rec_len
  CHARACTER(len=80) :: header(1000)
  CHARACTER(len=filenamelen) :: eigenvec_file, eigenval_file, maskfile=''
  CHARACTER(len=filenamelen) :: nlbias_file, argument
  LOGICAL :: there, noInvert=.FALSE.
  REAL(dp), POINTER :: plm(:, :)


  IF (nArguments() < 5)  THEN
     WRITE (*,'(/,a,/)') &
          'Usage: evecs2bias <eigenval_file> <eigenvec_file> <nside> '//&
          '<nmaps> <nlbias_file> [--noinv] [--limit <limit>] '//&
          '[-mask <maskfile>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     eigenval_file = TRIM(argument)
     INQUIRE(file=TRIM(eigenval_file), exist=there)
     IF (.NOT. there) THEN
        WRITE(*,*) TRIM(eigenval_file) // ' does not exist!'
        STOP
     END IF
     WRITE (*,*) ' Reading ' // TRIM(eigenval_file)

     CALL getArgument(2, argument)
     eigenvec_file = TRIM(argument)
     INQUIRE(file=TRIM(eigenvec_file), exist=there)
     IF (.NOT. there) then
        WRITE (*,*) TRIM(eigenvec_file) // ' does not exist!'
        STOP
     END IF
     WRITE (*,*) ' Reading ' // TRIM(eigenvec_file)

     CALL getArgument(3, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside == ', nside

     CALL getArgument(4, argument)
     READ (argument, *) nmaps
     WRITE (*,'(a,i4)') 'Nmaps    == ', nmaps

     CALL getArgument(5, argument)
     nlbias_file = TRIM(argument)
     WRITE (*,*) ' Writing to ' // TRIM(nlbias_file)

     
     IF (nArguments() > 5) THEN
        iArgument = 6
        do
           CALL getArgument(iargument, argument)
           IF (INDEX(argument, '-no') /= 0) THEN
              noInvert = .TRUE.
              WRITE (*,*) 'Using the eigenvalues without inverse'
           else if (index(argument, '-mask') /= 0) then
              iArgument = iArgument + 1
              CALL getArgument(iArgument, argument)
              maskfile = TRIM(ADJUSTL(argument))
              inquire(file=trim(maskfile), exist=there)
              if (.not. there) stop 'maskfile does not exist'
              write (*,*) 'Applying mask from ',trim(maskfile)
           ELSE IF (INDEX(argument, '-lim') /= 0) THEN
              iArgument = iArgument + 1           
              CALL getArgument(iargument, argument)
              READ (argument,*) limit
              WRITE (*,'("Excluding eigenvalues below ",g20.10)') limit
           ELSE
              WRITE (*,*) 'Unrecognized argument: '//TRIM(argument)
              STOP
           END IF
           iArgument = iArgument + 1
           if (iArgument > nArguments()) exit
        end do
        iArgument = iArgument + 1
     END IF
  END IF


  CALL cpu_time(t0)

  Npix  = 12*nside**2
  NNpix = nmaps*Npix
  lmax  = 4*nside

  zbounds = 0

  PRINT*,''
  PRINT*,'lmin, lmax, nside, Npix, nmaps = ',0,lmax,nside,Npix,nmaps

  ALLOCATE(w8ring_TQU(1:2*nside, 1:nmaps))
  w8ring_TQU = 1

  ALLOCATE(eigenvalue(0:NNpix-1))
  ALLOCATE(dmapout(0:Npix-1,1:nmaps))

  if (len_trim(maskfile) > 0) call input_mask(maskfile, mask, npix, .true.)

  ncl = 1
  if (nmaps == 3) ncl = 6
  ALLOCATE(nl(0:lmax, 1:ncl))
  ALLOCATE(alm_TGC_test(1:nmaps, 0:lmax, 0:lmax))
  ALLOCATE(nltest(0:lmax, 1:ncl))

  ALLOCATE(eigenvector(0:NNpix-1))

  nPlm = nside*(lMax+1)*(2*lMax-lMax+2)
  ALLOCATE(plm(0:nPlm-1, nMaps), stat=ierr)
  if (ierr /= 0) stop 'no room to precompute p_lm'
  CALL plm_gen(nSide, lMax, lMax, plm)

  PRINT*,''

  PRINT*,'Reading eigenvalues from  ',TRIM(ADJUSTL(eigenval_file))
  PRINT*,'Reading eigenvectors from ',TRIM(ADJUSTL(eigenvec_file))

  OPEN(3,file=eigenval_file)
  READ(3,*) eigenvalue(0:NNpix-1)
  CLOSE(3)

  if (limit == 0 .and. eigenvalue(0) < 0) then
     limit = abs(eigenvalue(0)/eigenvalue(NNpix-1))
     write (*,'("Found negative eigenvalue, set limit to ",es15.8)') limit
  end if
  where (eigenvalue <= limit*maxval(eigenvalue)) eigenvalue = 0

  IF (.NOT. noInvert) THEN
     write (*,*) 'Inverting the eigenvalues'
     where (eigenvalue /= 0) eigenvalue = 1.0/eigenvalue
  END IF

  PRINT*,'eigenvalue(0:2) = ',eigenvalue(0:2)
  PRINT*,'eigenvalue(NNpix-1) = ',eigenvalue(NNpix-1)

  PRINT*,''
  PRINT*,'Compute the noise bias'
  CALL cpu_time(t3)
  nl = 0.0d0

  inquire(iolength=rec_len) eigenvector
  open(unit=evec_unit, file=trim(eigenvec_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  DO i = 0,NNpix-1
     IF (eigenvalue(i) == 0) then
        write(*,*) 'Skipping mode #',i
        cycle
     end IF

     read (evec_unit, rec=i+1) eigenvector

     ! Subtract T mean
     eigenvector(:npix) = eigenvector(:npix) - sum(eigenvector(:npix)) / npix

     if (len_trim(maskfile) > 0) then
        do imap=0,nmaps-1
           eigenvector(imap*npix:(imap+1)*npix-1) = &
                eigenvector(imap*npix:(imap+1)*npix-1) * mask(:, 1)
        end do
     end if

     !if (any(isnan(eigenvector))) stop 'Eigenvector contains a NaN'
     
     dmapout = reshape(eigenvector, (/Npix, nmaps/))
     CALL convert_nest2ring(nside, dmapout)
     if (nmaps == 1) then
        CALL map2alm(nside,lmax,lmax,dmapout(:,1),alm_TGC_test,zbounds,&
             w8ring_TQU)
     else
        CALL map2alm(nside,lmax,lmax,dmapout,alm_TGC_test,zbounds,&
             w8ring_TQU,plm)
     end if
     CALL alm2cl(lmax,lmax,alm_TGC_test,nltest)

     !if (any(isnan(nltest))) then!stop 'Spectrum of the eigenvector contains a NaN'
     !   write (*,'("Warning: NaN in the spectrum for eigenmode ",i0)') i
     !   cycle
     !end if

     nl = nl + nltest*eigenvalue(i)

     IF (MODULO(i, NNpix/10)==0) WRITE (*,'(i4,a)') i*100/NNpix,' % done'
  END DO
  close(evec_unit)

  CALL cpu_time(t4)
  PRINT*,'Time to compute the noise_bias [s] = ',t4-t3

  !if (any(isnan(nl))) write (*,'("WARNING: Found NaN in the bias spectrum.")')

  OPEN(3,file=nlbias_file)
  DO l = 0,lmax
     if (ncl == 6) then
        WRITE(3,'(1I5,6E17.9)') l,nl(l,1),nl(l,2),nl(l,3),nl(l,4),nl(l,5),nl(l,6)
     else
        WRITE(3,'(1I5,6E17.9)') l,nl(l,1)
     end if
  END DO
  CLOSE(3)

  ! Write the noise bias into fits file
  header = ''
  call add_card(header, 'creator', 'evecs2bias', 'software creating the cls')
  call add_card(header, 'nmaps', nmaps, 'number of stokes parameters')
  call add_card(header, 'lmax', lmax, 'maximum multipole for cls')
  call add_card(header, 'units', 'K^2(ant)', 'unit of the cls')
  call add_card(header,"HISTORY","Input evals = "//TRIM(eigenval_file))
  call add_card(header,"HISTORY","Input evecs = "//TRIM(eigenvec_file))

  call add_card(header, 'TTYPE1', 'TT', 'label for field   1')
  if (ncl == 6) then
     call add_card(header, 'TTYPE2', 'EE', 'label for field   2')
     call add_card(header, 'TTYPE3', 'BB', 'label for field   3')
     call add_card(header, 'TTYPE4', 'TE', 'label for field   4')
     call add_card(header, 'TTYPE5', 'TB', 'label for field   5')
     call add_card(header, 'TTYPE6', 'EB', 'label for field   6')
  end if
     
  nlheader = SIZE(header)
  nlbias_file = '!' // trim(nlbias_file)
  if (index(nlbias_file,'.dat') > 0 .or. index(nlbias_file, '.fits') > 0) then
     nlbias_file(len_trim(nlbias_file)-3:len_trim(nlbias_file)+1) = '.fits'
  else
     nlbias_file = TRIM(nlbias_file) // '.fits'
  end if
  call write_asctab (nl, lmax, ncl, header, nlheader, nlbias_file)

contains

  subroutine input_mask(maskfile, mask, npix, verbose)

    use udgrade_nr, only : udgrade_nest

    character(len=*) :: maskfile
    real(sp), allocatable :: mask(:,:)
    integer(i4b) :: npix
    logical(lgt) :: verbose

    integer(i4b) :: npix_mask, nmaps_mask, ordering_mask, ierr, nside_mask, nside
    real(sp), allocatable :: masktemp(:,:)
    logical(lgt) :: there

    npix_mask = int(getsize_fits(maskfile, nmaps_mask, ordering_mask), i4b)
    if (verbose) then
       inquire(file=maskfile, exist=there)
       if (.not. there) stop 'maskfile does not exist!'
       write (*,'(a,i0)') ' npix_mask == ', npix_mask
       write (*,'(a,i0)') ' nmaps_mask == ', nmaps_mask
       select case (ordering_mask)
       case (0)
          write (*,'(a)') ' ordering_mask == UNKNOWN'
       case (1)
          write (*,'(a)') ' ordering_mask == RING'
       case (2)
          write (*,'(a)') ' ordering_mask == NEST'
       end select
    end if

    nmaps_mask = 1 ! For now, just use the first component
    
    allocate(masktemp(0:npix_mask-1, 1:nmaps_mask), stat=ierr)
    if (ierr /= 0) stop 'No room for masktemp'

    allocate(mask(0:npix-1, nmaps_mask), stat=ierr)
    if (ierr /= 0) stop 'No room for mask'

    nside_mask = npix2nside(npix_mask)
    call input_map(maskfile, masktemp, npix_mask, nmaps_mask)
    if (ordering_mask == 1) call convert_ring2nest(nside_mask, masktemp)
    
    nside = npix2nside(npix)
    call udgrade_nest(masktemp, nside_mask, mask, nside)
    
    deallocate(masktemp)
  end subroutine input_mask

END PROGRAM evecs2bias


