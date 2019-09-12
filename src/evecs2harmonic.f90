! This MPI program reads in the eigenvectors of a
!  pixel-pixel covariance and transforms them into a
! harmonic covariance matrix, where each element is a
! covariance of two harmonic coefficients: <a*_lm a_l'm'>
!
! the transform is computed by running map2alm
! on the columns and row of the matrix
!
! 2008-08-14 Reijo Keskitalo


PROGRAM covmat2harmonic

  USE healpix_types
  USE alm_tools
  USE pix_tools
  USE extension, ONLY : getArgument, nArguments

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(len=filenamelen) :: argument, evec_file, eval_file, outfile
  INTEGER, PARAMETER         :: covmat_unit=55
  REAL(dp), POINTER          :: evecs1(:,:), evals(:)
  COMPLEX(dpc), POINTER      :: evecs2(:,:)
  REAL(dp), POINTER          :: map(:,:)
  COMPLEX(dpc), POINTER      :: alm(:,:,:)
  LOGICAL                    :: there
  INTEGER(i4b)               :: ierr, irow, icol, id, ntasks, lmax, mmax, nalm
  INTEGER(i4b)               :: firstCol, lastCol, npix, nstokes, nside, isig
  INTEGER(i4b)               :: i, ell, m, firstIndex, lastIndex, itask, lmin
  INTEGER(i4b)               :: ncols_proc, nrows_proc, iroot, nalm2, count
  INTEGER(i4b)               :: nplm, nskip, ialm
  REAL(dp)                   :: zbounds(2), rcond_limit=1E-4
  REAL(dp), POINTER          :: w8ring(:,:), plm(:,:)
  COMPLEX(dpc), POINTER      :: sendbuf(:), recvbuf(:)

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (nArguments() < 6)  THEN
     IF (id == 0) WRITE(*,*) &
          'Usage: evecs2harmonic  '// &
          '<eigenvecs> <eigenvals> <nside> <nstokes> <lmin> <lmax> [<out>]'
     CALL mpi_finalize(ierr)
     STOP
  ELSE IF (id == 0) THEN
     CALL getArgument(1, argument)
     evec_file = TRIM(argument)
     INQUIRE(file=TRIM(evec_file), exist=there)
     IF (.NOT. there) STOP 'evec_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(evec_file)

     CALL getArgument(2, argument)
     eval_file = TRIM(argument)
     INQUIRE(file=TRIM(eval_file), exist=there)
     IF (.NOT. there) STOP 'eval_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(eval_file)

     CALL getArgument(3, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(4, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     CALL getArgument(5, argument)
     READ (argument, *) lmin
     WRITE (*,'(a,i4)') 'lMin    == ', lmin

     CALL getArgument(6, argument)
     READ (argument, *) lmax
     WRITE (*,'(a,i4)') 'lMax    == ', lmax

     IF (nArguments() > 6) THEN
        CALL getArgument(7, argument)
        outfile = TRIM(argument)
        WRITE (*,*) ' writing to ' // TRIM(outfile)
     END IF
  END IF
  
  ! Root node reads in the eigenvalues and determines how many eigenmodes
  ! to skip. 
  ! eigenvalues are assumed to correspond to the inverse covariance matrix
  CALL mpi_bcast(nside,    1,          MPI_INTEGER,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(evec_file,filenamelen,MPI_CHARACTER,0,mpi_comm_world,ierr)
  CALL mpi_bcast(nstokes,  1,          MPI_INTEGER,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(lmin,     1,          MPI_INTEGER,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(lmax,     1,          MPI_INTEGER,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(outfile,  filenamelen,MPI_CHARACTER,0,mpi_comm_world,ierr)
  npix = 12*nside**2
  allocate(evals(0:npix*nstokes-1), stat=ierr)
  if (ierr /= 0) stop 'no room for evals'
  if (id == 0) then
     open(unit=covmat_unit, file=trim(adjustl(eval_file)), status='old', &
          form='formatted')
     read(covmat_unit,*) evals
     close(covmat_unit)
     nskip = 0
     do
        if (evals(nskip)/evals(npix*nstokes-1) > rcond_limit) exit
        nskip = nskip + 1
     end do
     if (nskip > 0) write (*,*) ' Skipping ', nskip, 'bad eigenmodes'
     evals = 1.0/evals ! invert the eigenvalues for recomposition
  end if
  CALL mpi_bcast(evals,npix*nstokes,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
  CALL mpi_bcast(nskip, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)

  ! compute the number of nondegenerate alm coeffients
  nalm  = (lmax+1)**2 - lmin**2 ! all alms
  nalm2 = nalm - (nalm - (lmax+1-lmin))/2 ! exclude degenerate modes
  IF (id == 0) THEN
     WRITE (*,'(a,i6)') 'npix    == ', npix
     WRITE (*,'(a,i6)') 'nalm(lmin,lmax)     == ', nalm
     WRITE (*,'(a,i6)') 'nalm(nondegenerate) == ', nalm2
  END IF
  IF (MODULO(npix*nstokes,ntasks) /= 0 ) THEN
     IF (id == 0) WRITE (*,'(a,/,a)') &
          'Error: parallelization requires npix*nstokes ',&
          'to be divisible by number of processors'
     CALL mpi_finalize(mpi_comm_world, ierr)
     STOP
  END IF
  mmax = lmax

  ncols_proc = npix*nstokes/ntasks
  firstCol   = ncols_proc*id
  lastCol    = ncols_proc*(id+1) - 1

  ! read in the eigenvectors into local arrays
  ALLOCATE(evecs1(0:npix*nstokes-1, firstCol:lastCol), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat1'
  DO itask = 0, ntasks - 1
     IF (id == itask) THEN
        OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(evec_file)), status='old',&
             form='unformatted', access='direct', recl=npix*nstokes*8)
        DO icol = firstCol, lastCol
           READ(covmat_unit, rec=icol+1) evecs1(:, icol)
        END DO
        CLOSE(covmat_unit)
     END IF
     CALL mpi_barrier(mpi_comm_world, ierr)
  END DO

  ! convert the eigenvectors
  nplm = nside*(mmax + 1)*(2*lmax - mmax + 2)
  ALLOCATE(map(0:npix-1, nstokes), alm(nstokes, 0:lmax, 0:mmax), &
       w8ring(2*nside, nstokes), plm(0:nplm-1, nstokes), &
       evecs2(nstokes*nalm2, firstCol:lastCol), &
       stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat2'
  CALL plm_gen(nside, lmax, mmax, plm) ! Precompute the recursion factors
  zbounds = 0.0_dp
  w8ring  = 1.0_dp
  DO icol = firstCol, lastCol
     DO isig = 1, nstokes
        map(:, isig) = evecs1(npix*(isig-1):npix*isig-1, icol)
     END DO
     CALL convert_nest2ring(nside, map) ! covariance matrix is NESTED
     CALL map2alm(nside, lmax, mmax, map, alm, zbounds, w8ring, plm)
     count = 1
     DO isig = 1, nstokes
        DO ell = lmin, lmax
           DO m = 0, ell
              evecs2(count, icol) = alm(isig, ell, m)
              count = count + 1
           END DO
        END DO
     END DO
  END DO


  ! Then compose the harmonic covariance matrix from its
  ! eigenvectors one column at a time. Root node writes it
  ! to disk
  IF (id == 0) then
     WRITE (*,'(/,a,i4,a,i8,a)') &
       ' Saving the ', 3*nalm2,'^2 ==', (3*nalm2)**2, ' complex matrix'
     OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(outfile)), status='replace', &
             form='unformatted', access='direct', recl=2*8*nalm2*nstokes)
  END IF
  allocate(sendbuf(nalm2*nstokes), recvbuf(nalm2*nstokes), stat=ierr)
  if (ierr /= 0) stop 'No room for buffers'
  DO ialm = 1, nalm2*nstokes
     if (id == 0) recvbuf = 0
     sendbuf = 0
     ! compute the same column of \vec u_i \vec u_i^\dagger
     ! for all eigenvectors and sum the results
     do icol = firstCol, lastCol
        if (icol < nskip) cycle
        sendbuf = sendbuf + evals(icol)*evecs2(:,icol)*conjg(evecs2(ialm,icol))
     end do

     call mpi_reduce(sendbuf, recvbuf, nalm2*nstokes, MPI_DOUBLE_COMPLEX, &
          MPI_SUM, 0, mpi_comm_world, ierr)

     if (id == 0) write (covmat_unit, rec=ialm) recvbuf
  END DO
  if (id == 0) CLOSE(covmat_unit)

  CALL mpi_finalize(ierr)


END PROGRAM covmat2harmonic
