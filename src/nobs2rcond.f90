! Compute the reciprocal condition number of the supplied 3x3 matrices and save
! it as a map.
!
!  04-08-2010 RK
!
!  Changes:
!     2013-06-12 : Buffered processing to reduce memory requirement
!                   and overflow of 32 bit indices in Healpix

PROGRAM nobs2rcond

  USE healpix_types
  USE fitstools, ONLY   : output_map, getsize_fits, input_map, input_tod
  USE head_fits, ONLY   : add_card, write_minimal_header
  USE pix_tools, ONLY   : convert_ring2nest, nside2npix
  USE rngmod
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments
  use omp_lib

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument
  CHARACTER(len=filenamelen) :: file_nobs, file_map
  CHARACTER(len=80) :: header(1000)
  INTEGER(i4b) :: nmaps, ncols, nside
  integer(i8b) :: npix_tot, npix, ipix
  INTEGER(i4b) :: iargument, ordering, info
  INTEGER, PARAMETER :: covmat_unit=55, buflen=1e5
  REAL(sp), POINTER :: rcond(:, :)
  REAL(dp), POINTER :: nobs(:, :)
  REAL(dp), POINTER :: eigenvals(:), eigenvecs(:, :)
  LOGICAL :: invert_matrix=.false.
  INTEGER(i4b) :: ierr
  integer(i8b) :: workspacelen=100, firstpix, nread, ipixtot
  real(dp), pointer :: workspace(:)
  INTEGER(i4b) :: row, col, i

  IF (nArguments() == 0)  THEN
     write (*,'(/,a,/)') ' Usage: nobs2rcond <file_nobs> [--inv] [-o <file_map>]'
     stop
  ELSE
     CALL getArgument(1, file_nobs)
     WRITE (*,*) ' file_nobs = '//trim(file_nobs)

     file_map = '!noisemap.fits'
     invert_matrix = .false.
     iargument = 2
     do
        if (nArguments() < iargument) exit
        CALL getArgument(iargument, argument)

        if (index(argument, '-o') > 0) then
           iargument = iargument + 1
           call getArgument(iargument, file_map)
        else if (index(argument, '-inv') > 0) then
           invert_matrix = .true.
           write (*,*) 'Inverting the nobs matrices'
        else
           write (*,*) 'ERROR: unrecognized option: '//trim(argument)
           stop
        end if
        iargument = iargument + 1
     end do

     WRITE (*,*) ' Generating ' // TRIM(file_map)
  END IF

  npix_tot = getsize_fits(file_nobs, nmaps=ncols, ordering=ordering, nside=nside)
  !nmaps = floor( sqrt(ncols * 2.0) )
  nmaps = nint(0.5 * ( sqrt(1.0 + 8 * ncols ) - 1 ))
  npix = nside2npix(nside)

  write (*,'(" nmaps == ",i8)') nmaps
  write (*,'(" nside == ",i8)') nside
  write (*,'("  npix == ",i8)') npix

  ALLOCATE(rcond(0:npix-1,1), nobs(buflen,ncols), stat=ierr)
  IF (ierr /= 0) STOP 'No room for map or nobs'

  !CALL tic
  !call input_map(file_nobs, nobs, npix, ncols)
  !call input_tod(file_nobs, nobs, int(npix,i8b), ncols, firstpix=0_i8b)
  !call toc('Read nobs')

  call tic

  do firstpix = 0,npix,buflen

     nread = buflen
     if (firstpix + nread > npix) nread = npix - firstpix

     call input_tod(file_nobs, nobs(1:nread,:), nread, ncols, firstpix=firstpix)

     !$OMP PARALLEL DEFAULT(PRIVATE) &
     !$OMP    SHARED(npix, nmaps, nobs, rcond, workspacelen, nread, firstpix)

     allocate(workspace(workspacelen), eigenvecs(nmaps,nmaps), &
          eigenvals(nmaps), stat=ierr)

     !$OMP DO SCHEDULE(STATIC,4)
     DO ipix = 1, nread
        ipixtot = ipix + firstpix - 1
        if (mod(ipixtot,npix/10) == 0) then
           write (*,'(1x,i3,"%")') nint(ipixtot*10._dp / npix)*10
        end if

        if (nobs(ipix, 1) <= 0) then
           rcond(ipixtot, 1) = HPX_SBADVAL
           cycle
        end if

        eigenvecs(1, 1) = nobs(ipix, 1)

        if (nmaps > 1) then
           i = 0
           do row = 1, nmaps
              do col = row, nmaps
                 i = i + 1
                 eigenvecs(row, col) = nobs(ipix, i)
                 if (row /= col) eigenvecs(col, row) = nobs(ipix, i)
              end do
           end do

           ! compute the eigenvalue decomposition
           ! http://www.netlib.org/lapack/double/dsyev.f
           call dsyev('V', 'U', nmaps, eigenvecs, nmaps, eigenvals, workspace, workspacelen, info)
           if (info == 0) then
              rcond(ipixtot, 1) = real(minval(eigenvals)/maxval(eigenvals), sp)
           else
              if (info < 0) write (*,'("dsyev failed. info == ",i0)') info
              rcond(ipixtot, 1) = HPX_SBADVAL
           end if
        else
           rcond(ipixtot, 1) = 1
        end if
     END DO
     !$OMP END DO
     deallocate(workspace, eigenvecs, eigenvals)
     !$OMP END PARALLEL
  end do

  CALL toc('Generate map')

  CALL write_minimal_header(header, 'MAP', nside=nside, ordering='NESTED', &
       creator='nobs2rcond', polar=.false.)

  CALL output_map(rcond, header, file_map)


END PROGRAM nobs2rcond
