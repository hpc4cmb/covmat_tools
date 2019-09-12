! This code reads in the eigenvectors and eigenvalues
! and composes the matrix they correspond to
!
! Revisions:
!   2012-10-25 : Fixed a serious bug when handling downgraded (non-square)
!                eigenvector inputs
 
PROGRAM evecs2covmat

  USE healpix_types
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY : getArgument, nArguments
  use fitstools, only : input_map
  use pix_tools, only : nside2npix

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, evec_file, eval_file, outfilename
  INTEGER(i4b) :: nstokes, nside_out, nside_in
  INTEGER(i8b) :: icol, npix_in, npix_out
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:, :), evec(:, :), evecT(:, :), eval(:)
  REAL(dp), POINTER :: buffer(:)
  REAL(dp) :: evalmax, limit=1E-30
  LOGICAL(lgt) :: there, invert=.false., sqroot=.false.
  INTEGER(i4b) :: ierr, iArgument, nskip, nthreads
  INTEGER(i4b) :: omp_get_num_threads, rec_len

  nthreads = omp_get_num_threads()
  write (*,*) 'Running with ',nthreads,' threads'

  IF (nArguments() < 4)  THEN
     WRITE (*,'(/,a,/)') 'Usage: evecs2covmat <evals> <evecs> <nside_in> ' // &
          '<nside_out> <nstokes> [--inv] [--sqrt] [-o <outfile>] [-lim <limit>] '
     STOP
  ELSE
     CALL getArgument(1, argument)
     eval_file = TRIM(argument)
     INQUIRE(file=TRIM(eval_file), exist=there)
     IF (.NOT. there) STOP 'eval_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(eval_file)

     CALL getArgument(2, argument)
     evec_file = TRIM(argument)
     INQUIRE(file=TRIM(evec_file), exist=there)
     IF (.NOT. there) STOP 'evec_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(evec_file)

     CALL getArgument(3, argument)
     READ (argument, *) nside_in
     WRITE (*,'(a,i4)') 'Nside_in   == ', nside_in

     CALL getArgument(4, argument)
     READ (argument, *) nside_out
     WRITE (*,'(a,i4)') 'Nside_out  == ', nside_out

     CALL getArgument(5, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     outfilename = 'composed_covmat.dat'
     iArgument   = 6
     do
        CALL getArgument(iArgument, argument)
        if (index(argument, '-o') /= 0) then
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           outfilename = TRIM(ADJUSTL(argument))
        else if (index(argument, '-lim') /= 0) then
           iArgument = iArgument + 1
           CALL getArgument(iArgument, argument)
           read(argument, *, iostat=ierr) limit
           if ( ierr /= 0 ) then
              print *,'Failed to parse limit from ' // trim(argument)
              stop
           end if
           write (*,*) 'Using threshold eval > evalmax * ',limit
        else if (index(argument, '-inv') /= 0) then
           invert = .true.
           write (*,*) 'Inverting while composing'
        else if (index(argument, '-sqrt') /= 0) then
           sqroot = .true.
           write (*,*) 'Composing the square root'
        else
           write (*,*) 'Unrecognized command line option: ',trim(argument)
           stop
        end if

        iArgument = iArgument + 1
        if (iArgument > nArguments()) exit
     end do

     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix_in = nside2npix(nside_in)
  npix_out = nside2npix(nside_out)
  write(*,'(" npix_in == ",i0)') npix_in
  write(*,'(" npix_out == ",i0)') npix_out

  ! read, multiply and write. one column at a time
  CALL tic
  ALLOCATE(covmat(0:npix_out*nstokes-1,0:npix_out*nstokes-1), &
       eval(0:npix_in*nstokes-1), evec(0:npix_out*nstokes-1,0:npix_in*nstokes-1), &
       evecT(0:npix_in*nstokes-1,0:npix_out*nstokes-1), buffer(npix_out*nstokes), &
       stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'

  OPEN(unit=covmat_unit, file=TRIM(eval_file), status='old', form='formatted')
  read(covmat_unit,*) eval
  CLOSE(covmat_unit)

  ! process the eigenvalues
  evalmax = maxval(eval)
  limit = evalmax*limit
  where(eval < limit) eval = 0.0_dp
  if (sqroot) where(eval > 0.0_dp) eval = sqrt(eval)
  if (invert) where(eval > 0.0_dp) eval = 1.0_dp/eval
  nskip = count(eval == 0.0_dp)

  inquire(iolength=rec_len) buffer
  OPEN(unit=covmat_unit+1, file=TRIM(evec_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  do icol=0,npix_in*nstokes-1
     read(covmat_unit+1, rec=icol+1) buffer
     evec(:,icol) = buffer
  end do
  CLOSE(covmat_unit+1)

  covmat = 0.0
  eval = sqrt(eval)
  DO icol = 0, npix_in*nstokes - 1
     evec(:,icol) = evec(:,icol)*eval(icol)
  END DO
  write (*,'(a,i0,a,i0)') ' Omitted ', nskip, ' modes, kept ', &
       npix_in*nstokes-nskip

  evecT = transpose(evec)
  covmat = matmul(evec, evecT)
  ! SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
  !   $                   BETA, C, LDC )
  ! C := alpha*op( A )*op( B ) + beta*C,
  !n = npix*nstokes
  !call sgemm('N','T', n, n, n, 1, evec, n, evec, n, 0, covmat, n)

  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit+2, file=TRIM(outfilename), status='replace', &
       form='unformatted', access='direct', recl=rec_len)
  write(covmat_unit+2, rec=1) covmat
  CLOSE(covmat_unit+2)

  CALL toc('compose_covmat')

  deallocate(covmat, eval, evec, evecT, buffer)

END PROGRAM evecs2covmat
