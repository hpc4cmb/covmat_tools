! This MPI program reads in a binary pixel-pixel covariance
! matrix in block (II, IQ, ..) format and transforms it into
! harmonic space, where each element is a covariance of two
! harmonic coefficients: <a*_lm a_l'm'>
!
! the transform, Y* N Y is computed by running map2alm
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

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfile
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat1(:, :), covmat2(:, :), covmat3(:, :)
  COMPLEX(dpc), POINTER :: covmat4(:, :)
  REAL(dp), POINTER :: map(:, :)
  COMPLEX(dpc), POINTER :: alm(:, :, :)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, irow, icol, id, ntasks, lmax, mmax, nalm
  INTEGER(i4b) :: firstCol, lastCol, npix, nstokes, nside, isig
  INTEGER(i4b) :: ell, m, firstIndex, lastIndex, itask, lmin
  INTEGER(i4b) :: ncols_proc, nrows_proc, iroot, nalm2, count
  INTEGER(i4b) :: nplm, rec_len
  REAL(dp) :: zbounds(2)
  REAL(dp), POINTER :: w8ring(:, :), sendbuf(:), recvbuf(:), plm(:, :)

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (nArguments() < 5)  THEN
     IF (id == 0) WRITE(*,*) &
          'Usage: covmat2harmonic  '// &
          '<covmat> <nside> <nstokes> <lmin> <lmax> [<out>]'
     CALL mpi_finalize(ierr)
     STOP
  ELSE IF (id == 0) THEN
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     CALL getArgument(4, argument)
     READ (argument, *) lmin
     WRITE (*,'(a,i4)') 'lMin    == ', lmin

     CALL getArgument(5, argument)
     READ (argument, *) lmax
     WRITE (*,'(a,i4)') 'lMax    == ', lmax

     IF (nArguments() > 5) THEN
        CALL getArgument(6, argument)
        outfile = TRIM(argument)
        WRITE (*,*) ' writing to ' // TRIM(outfile)
     END IF
  END IF

  CALL mpi_bcast(covmat_file,filenamelen,mpi_character,0,mpi_comm_world,ierr)
  CALL mpi_bcast(nside,      1,          mpi_integer,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(nstokes,    1,          mpi_integer,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(lmin,       1,          mpi_integer,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(lmax,       1,          mpi_integer,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(outfile,    filenamelen,mpi_character,0,mpi_comm_world,ierr)

  npix  = 12*nside**2
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
  firstCol = ncols_proc*id
  lastCol  = ncols_proc*(id+1) - 1

  ! read in the covariance matrix
  ALLOCATE(covmat1(0:npix*nstokes-1, firstCol:lastCol), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat1'
  DO itask = 0, ntasks - 1
     IF (id == itask) THEN
        inquire(iolength=rec_len) covmat1(:,firstCol)
        OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(covmat_file)), status='old',&
             form='unformatted', access='direct', recl=rec_len)
        DO icol = firstCol, lastCol
           READ(covmat_unit, rec=icol+1) covmat1(:, icol)
        END DO
        CLOSE(covmat_unit)
     END IF
     CALL mpi_barrier(mpi_comm_world, ierr)
  END DO

  ! convert the columns of the matrix, this is the (Y N) part
  nplm = nside*(mmax + 1)*(2*lmax - mmax + 2)
  ALLOCATE(map(0:npix-1, nstokes), alm(nstokes, 0:lmax, 0:mmax), &
       w8ring(2*nside, nstokes), plm(0:nplm-1, nstokes), &
       covmat2(nstokes*nalm2*2, firstCol:lastCol), &
       stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat2'
  CALL plm_gen(nside, lmax, mmax, plm) ! Precompute the recursion factors
  zbounds = 0.0_dp
  w8ring  = 1.0_dp
  DO icol = firstCol, lastCol
     DO isig = 1, nstokes
        map(:, isig) = covmat1(npix*(isig-1):npix*isig-1, icol)
     END DO
     CALL convert_nest2ring(nside, map)
     CALL map2alm(nside, lmax, mmax, map, alm, zbounds, w8ring, plm)
     ! save the required elements of the complex vector into
     ! a real matrix. Every second row will be the complex part
     ! of the previous row
     count = 1
     DO isig = 1, nstokes
        DO ell = lmin, lmax
           DO m = 0, ell
              ! complex conjugate the alms, since map2alm returns
              ! sum_p Y*_lm(p) map(p) = a_lm
              covmat2(count,     icol) = REAL(alm(isig, ell, m), dp)
              covmat2(count + 1, icol) = -AIMAG(alm(isig, ell, m))
              count = count + 2
           END DO
        END DO
     END DO
  END DO


  ! transpose the matrix, to make this part smooth we need
  ! npix*nstokes to be divisible by ntasks
  nrows_proc = 2*INT(nalm2*nstokes/ntasks, i4b)
  firstIndex = nrows_proc*id + 1
  lastIndex  = nrows_proc*(id+1)
  if (id == ntasks-1) lastIndex = nalm2*nstokes*2
  ALLOCATE(covmat3(0:npix*nstokes-1, firstIndex:lastIndex), &
       sendbuf(ncols_proc), recvbuf(npix*nstokes), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat3'
  DO irow = 1, nstokes*nalm2*2
     iroot   = min((irow-1)/nrows_proc, ntasks-1)
     sendbuf = covmat2(irow, :)
     CALL mpi_gather(sendbuf, ncols_proc, MPI_DOUBLE_PRECISION, &
          recvbuf, ncols_proc, MPI_DOUBLE_PRECISION, &
          iroot, mpi_comm_world, ierr)
     IF (id == iroot) covmat3(:, irow) = recvbuf
  END DO


  ! convert the columns of the matrix again to get Y N Y*^T
  ! Because of the transpose we end up having the transpose of the 
  ! alm-covariance matrix which we know to be hermitian. This is fixed by 
  ! taking the complex conjugate of the matrix
  ALLOCATE(covmat4(nstokes*nalm2, firstIndex/2+1:lastIndex/2), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat4'
  DO icol = firstIndex, lastIndex
     DO isig = 1, nstokes
        map(:, isig) = covmat3(npix*(isig-1):npix*isig-1, icol)
     END DO
     CALL convert_nest2ring(nside, map)
     CALL map2alm(nside, lmax, mmax, map, alm, zbounds, w8ring, plm)
     count = 1
     DO isig = 1, nstokes
        DO ell = lmin, lmax
           DO m = 0, ell
              IF (MODULO(icol,2) == 1) THEN
                 ! real part
                 covmat4(count, icol/2+1) = alm(isig, ell, m)
              ELSE
                 ! complex part to the previous column
                 covmat4(count, icol/2) = covmat4(count, icol/2) + CMPLX(&
                      -AIMAG(alm(isig,ell,m)), REAL(alm(isig, ell, m), dp), dpc)
              END IF
              count = count + 1
           END DO
        END DO
     END DO
  END DO

  ! complex conjugate the matrix to get N, not N^T=N*
  covmat4 = conjg(covmat4)


  ! save the matrix
  IF (id == 0) then
     WRITE (*,'(/,a,i4,a,i8,a)') &
       ' Saving the ', nalm2,'^2 ==', nalm2**2, ' complex matrix'
     inquire(iolength=rec_len) covmat4(:,firstIndex/2+1)
     OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(outfile)), status='replace', &
             form='unformatted', access='direct', recl=rec_len)
     CLOSE(covmat_unit)
  END IF
  ! Have each process write its columns one by one to the same common file.
  DO itask = 0, ntasks - 1
     IF (id == itask) THEN
        OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(outfile)), status='old', &
             form='unformatted', access='direct', recl=rec_len)
        DO icol = firstIndex/2+1, lastIndex/2
           WRITE (covmat_unit, rec=icol) covmat4(:, icol)
           WRITE (*,'("(",ES10.3,",",ES10.3,")")') covmat4(icol, icol)
        END DO
        CLOSE(covmat_unit)
     END IF
     call mpi_barrier(mpi_comm_world, ierr)
  END DO

  CALL mpi_finalize(ierr)


END PROGRAM covmat2harmonic
