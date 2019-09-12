! Generate white noise map according to the 3x3 nobs matrices
!
!  22-02-2010 RK

PROGRAM generate_whitenoisemap

  USE healpix_types
  USE fitstools, ONLY   : output_map, getsize_fits, input_map, input_tod
  USE head_fits, ONLY   : add_card, write_minimal_header
  USE pix_tools, ONLY   : convert_ring2nest, nside2npix
  USE rngmod
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, file_nobs, file_map
  CHARACTER(len=80) :: header(1000)
  INTEGER(i4b) :: nmaps, ncols, nside, isig
  integer(i8b) :: npix_tot, npix, ipix
  INTEGER(i4b) :: seed, iargument, ordering, info
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(sp), POINTER :: map(:, :)
  REAL(dp), POINTER :: nobs(:, :), covmat(:, :), iqu(:)
  REAL(dp), POINTER :: eigenvals(:), eigenvecs(:,:)
  LOGICAL :: there, invert_matrix=.false.,nobs_matrix = .false., pol=.false.
  INTEGER(i4b) :: ierr, ieigen
  TYPE(planck_rng) :: rnghandle
  REAL(dp) :: sigma_TT, sigma_pol
  integer(i8b), parameter :: workspacelen=100
  real(dp) :: workspace(workspacelen)

  IF (nArguments() < 2)  THEN
     write (*,'(/,a,/)') ' Usage: generate_whitenoisemap <nside> ' // &
          '[-sigma_TT <TT sigma> -sigma_pol <Pol sigma>] [-nobs <file_nobs>] ' // &
          '[--invert_matrices] [-o <file_map>]'
     stop
  ELSE
     CALL getArgument(1, argument)
     READ (argument, *, iostat=ierr) nside
     WRITE (*,*) 'nside == ', nside

     iargument = 2
     file_map = '!noisemap.fits'
     invert_matrix = .false.
     nobs_matrix = .false.
     pol=.false.
     DO
         IF (iargument > nArguments()) EXIT

         CALL getArgument(iargument,argument)
         IF (INDEX(argument, '-sigma_TT') /= 0) THEN
             iArgument = iArgument + 1 
             CALL getArgument(iargument, argument)
             READ (argument,*,iostat=ierr) sigma_TT
             IF (ierr /= 0) STOP 'Unable to parse sigma_TT'
             WRITE (*,*) 'sigma_TT == ', sigma_TT
         ELSE IF (INDEX(argument, '-sigma_pol') /= 0) THEN
             iArgument = iArgument + 1
             CALL getArgument(iargument, argument)
             pol = .true.
             READ (argument,*,iostat=ierr) sigma_pol
             IF (ierr /= 0) STOP 'Unable to parse sigma_pol'
             WRITE (*,*) 'sigma_pol == ', sigma_pol
         ELSE IF (INDEX(argument, '-nobs') /= 0) THEN
             iArgument = iArgument + 1
             CALL getArgument(iargument, argument)
             nobs_matrix = .true.
             file_nobs = TRIM(argument)
             INQUIRE(file=TRIM(file_nobs), exist=there)
             IF (.NOT. there) STOP 'file_nobs does not exist'
         ELSE IF (INDEX(argument, '-o') > 0) THEN
             iArgument = iArgument + 1
             CALL getArgument(iargument, file_map)
             INQUIRE(file=TRIM(file_map), exist=there)
             IF (there) STOP 'file_map exists'
         ELSE IF (INDEX(argument, '-inv') > 0) THEN
             invert_matrix = .true.
             write (*,*) 'Inverting the nobs matrices'
         ELSE
             write (*,*) 'ERROR: unrecognized option: '//trim(argument)
             stop
         END IF
         iargument = iargument + 1
     END DO

     WRITE (*,*) ' Generating ' // TRIM(file_map)
  END IF

  if (nobs_matrix) then
      npix_tot = getsize_fits(file_nobs, nmaps=ncols, ordering=ordering, &
           nside=nside)
      nmaps = 1; if (ncols == 6) nmaps = 3
  else if (pol) then
      nmaps = 3
  else
      nmaps = 1
  end if

  npix = nside2npix(nside)

  write (*,'(" nmaps == ",i8)') nmaps
  write (*,'(" nside == ",i8)') nside
  write (*,'("  npix == ",i8)') npix

  ALLOCATE(map(0:npix-1,nmaps), iqu(nmaps), stat=ierr)
  IF (ierr /= 0) STOP 'No room for map'

  IF (nobs_matrix) THEN
      ALLOCATE(nobs(0:npix-1,ncols), covmat(nmaps,nmaps), &
      eigenvals(nmaps), eigenvecs(nmaps,nmaps), stat=ierr)
      IF (ierr /= 0) STOP 'No room for map or nobs'

      CALL tic
      call input_map(file_nobs, nobs, npix, ncols)
      call toc('Read nobs')

  END IF

  call tic
  CALL SYSTEM_CLOCK(count=seed)
  CALL rand_init(rnghandle, seed)
  DO ipix = 0, npix-1

     IF (nobs_matrix) THEN
         if (nobs(ipix, 1) <= 0) then
            map(ipix,:) = HPX_SBADVAL
            cycle
         end if
     END IF

     ! generate the uniform gaussian iqu
     DO isig = 1, nmaps
        iqu(isig) = rand_gauss(rnghandle)
     END DO

     IF (nobs_matrix) THEN
         covmat(1,1) = nobs(ipix,1)
 
         if (nmaps == 3) then
            ! compute the cholesky decomposition of the 3x3 matrix
            covmat(1,2) = nobs(ipix,2)
            covmat(1,3) = nobs(ipix,3)
            covmat(2,1) = nobs(ipix,2)
            covmat(2,2) = nobs(ipix,4)
            covmat(2,3) = nobs(ipix,5)
            covmat(3,1) = nobs(ipix,3)
            covmat(3,2) = nobs(ipix,5)
            covmat(3,3) = nobs(ipix,6)
        
            if (invert_matrix) then
               ! compute the eigenvalue decomposition
               ! http://www.netlib.org/lapack/double/dsyev.f
               eigenvecs = covmat
               call dsyev('V', 'U', nmaps, eigenvecs, nmaps, eigenvals, &
                    workspace, workspacelen, info)
               if (info == 0) then
                  do ieigen = 1,nmaps
                     eigenvecs(:,ieigen) = &
                          eigenvecs(:,ieigen) / eigenvals(ieigen)**0.25 ! get the sqrt of the 3x3 matrix
                  end do
                  covmat = matmul(eigenvecs, transpose(eigenvecs))
                  iqu = matmul(covmat, iqu)
               else
                  if (info < 0) write (*,'("dsyev failed. info == ",i0)') info
                  iqu = 0
               end if
            else
               ! compute the cholesky decomposition
               ! http://www.netlib.org/lapack/lapack-3.1.1/SRC/dpotrf.f
               call dpotrf('L', nmaps, covmat, nmaps, info)
               if (info == 0) then
                   ! successful decomposition
                   covmat(1,2:3) = 0
                   covmat(2,3) = 0
                  iqu = matmul(covmat, iqu) ! CHECK THIS LINE
               else
                   if (info < 0) write (*,'("dpotrf failed. info == ",i0)') info
                   iqu = 0
               end if
            end if
         else
            if (invert_matrix) then
               if (covmat(1,1) /= 0) covmat(1,1) = 1.0/covmat(1,1)
            end if
            iqu(1) = iqu(1) * covmat(1,1)
         end if
     ELSE IF (pol) THEN
         iqu(1) = iqu(1) * sigma_TT
         iqu(2) = iqu(2) * sigma_pol 
         iqu(3) = iqu(3) * sigma_pol
     ELSE
         iqu(1) = iqu(1) * sigma_TT
     END IF
     map(ipix, :) = real(iqu, sp)
  END DO
  
  CALL toc('Generate map')

  CALL write_minimal_header(header, 'MAP', nside=nside, ordering='NESTED', &
       creator='generate_whitenoisemap', randseed=seed, polar=(nmaps>1))

  CALL output_map(map, header, file_map)


END PROGRAM generate_whitenoisemap
