
PROGRAM harmonic2bias

  USE healpix_types
  USE head_fits
  USE fitstools
  USE extension, ONLY   : nArguments, getArgument
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file
  CHARACTER(len=filenamelen) :: outfilename
  INTEGER(i4b) :: nstokes, lmin, lmax, nalm, nalm2, nell, ell, m
  REAL(dp) :: A
  INTEGER, PARAMETER :: covmat_unit=55
  COMPLEX(dpc), POINTER :: covmat(:)
  REAL(dp), POINTER :: bias(:,:)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, icol, nlheader, ncl
  CHARACTER(len=80) :: header(1000)


  IF (nArguments() < 4)  THEN
     STOP 'Usage: harmonic2angular <covmat> <lmin> <lmax> <nstokes> [<outfile>]'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) lmin
     WRITE (*,'(a,i4)') 'lmin    == ', lmin

     CALL getArgument(3, argument)
     READ (argument, *) lmax
     WRITE (*,'(a,i4)') 'lmax    == ', lmax

     CALL getArgument(4, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes  == ', nstokes

     IF (nArguments() == 5) THEN
        CALL getArgument(5, argument)
        outfilename = TRIM(ADJUSTL(argument))
     ELSE
        outfilename = 'bias.dat'
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  nalm  = (lmax+1)**2 - lmin**2 ! all alms
  nalm2 = nalm - (nalm - (lmax+1-lmin))/2 ! exclude degenerate modes
  nell = lmax - lmin + 1

  ncl = 1
  IF (nstokes == 3) ncl = 6
  ALLOCATE(covmat(nalm2*nstokes), bias(0:lmax, ncl), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat'

  CALL tic
  OPEN(unit=covmat_unit,   file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=nalm2*nstokes*2*8)
  bias = 0
  ell = lmin
  m   = 0
  DO icol = 1, nalm2*nstokes
     ! only the l=l' and m=m' elements count for th bias
     A = 2.0
     if (m == 0) A = 1.0
     READ (covmat_unit, rec=icol) covmat
     IF (icol <= nalm2) THEN
        bias(ell, 1) = bias(ell, 1) + A*real(covmat(icol))         ! TT bias
     ELSE IF (icol <= 2*nalm2) THEN
        bias(ell, 2) = bias(ell, 2) + A*real(covmat(icol))         ! EE
        bias(ell, 4) = bias(ell, 4) + A*real(covmat(icol-nalm2))   ! TE
     ELSE
        bias(ell, 3) = bias(ell, 3) + A*real(covmat(icol))         ! BB
        bias(ell, 5) = bias(ell, 5) + A*real(covmat(icol-nalm2*2)) ! TB
        bias(ell, 6) = bias(ell, 6) + A*real(covmat(icol-nalm2))   ! EB
     END IF
     m = m + 1
     if (m > ell) then
        ell = ell + 1
        if (ell > lmax) ell = lmin
        m = 0
     end if     
  END DO
  CLOSE(unit=covmat_unit)
  ! normalize the C_ells
  DO ell = 0, lmax
     bias(ell, :) = bias(ell, :) / (2*ell + 1)
  END DO
  CALL toc('convert covmat')

  CALL tic
  OPEN(covmat_unit+1,file=outfilename,status='replace',form='formatted')
  DO ell = 0, lmax
     IF (ncl == 6) THEN
        WRITE(covmat_unit+1, '(I5,6E17.9)') ell, bias(ell,:)
     ELSE
        WRITE(covmat_unit+1, '(I5,E17.9)') ell, bias(ell,:)
     END IF
  END DO
  CLOSE(covmat_unit+1)
  ! Write the noise bias into fits file
  header = ''
  CALL add_card(header, 'creator', 'harmonic2bias','software creating the cls')
  CALL add_card(header, 'nmaps', nstokes, 'number of stokes parameters')
  CALL add_card(header, 'lmin', lmin, 'minimum multipole for cls')
  CALL add_card(header, 'lmax', lmax, 'maximum multipole for cls')
  CALL add_card(header, 'units', 'K^2(ant)', 'unit of the cls')
  CALL add_card(header,"HISTORY","Input NCVM = "//TRIM(covmat_file))
  nlheader    = SIZE(header)
  outfilename = '!' // TRIM(outfilename)
  IF (INDEX(outfilename,'.dat') > 0) THEN
     outfilename(LEN_TRIM(outfilename)-3:LEN_TRIM(outfilename)+1) = '.fits'
  ELSE
     outfilename = TRIM(outfilename) // '.fits'
  END IF
  CALL write_asctab(bias, lmax, ncl, header, nlheader, outfilename)
  CALL toc('save C_ell covmat')

END PROGRAM harmonic2bias
