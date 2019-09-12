!> Small utility module for covariance matrices
!! 2007-09-07 RK
MODULE covmat_util

  USE healpix_types
  !USE scalapack_tools

  IMPLICIT NONE

  INTEGER, PARAMETER, private :: covmat_util_unit = 613 ! output unit
  REAL(dp), private :: time_start, tic_times(1000)
  INTEGER, private :: ID = 0

  INTERFACE eigen_invert_matrix
     MODULE PROCEDURE eigen_invert_symmetric_matrix, &
          eigen_invert_hermitian_matrix
  END INTERFACE



CONTAINS



  FUNCTION get_file_size_dp(file_name)

    ! return the size of the file in bytes

    use iso_c_binding

    IMPLICIT NONE

    EXTERNAL fsize_c

    CHARACTER(len=*) :: file_name
    INTEGER(i8b) :: get_file_size_dp

    character(c_char) :: c_name(1024)
    integer(c_int) :: c_name_len, i
    integer(c_long_long) :: c_size
    LOGICAL :: there
    
    INQUIRE(file=file_name, exist=there)
    IF (.NOT. there) THEN
       get_file_size_dp = -1
       RETURN
    END IF

    !call fsize_c(file_name, len_trim(file_name), get_file_size_dp)
    c_name_len = len_trim(file_name)
    do i = 1, c_name_len
       c_name(i) = file_name(i:i)
    end do
    call fsize_c(c_name, c_name_len, c_size)

    get_file_size_dp = c_size / 8_i8b ! convert to double elements

  END FUNCTION get_file_size_dp



  SUBROUTINE get_matrix_size_dp(file_name, nside, nstokes, npixtot, nelem)

    ! gets the file size and converts into healpix parameters

    IMPLICIT NONE

    CHARACTER(len=filenamelen), intent(in) :: file_name
    INTEGER(i4b), intent(out) :: nside, nstokes
    INTEGER(i8b), OPTIONAL, intent(out) :: npixtot
    INTEGER(i8b), OPTIONAL, intent(out) :: nelem

    INTEGER(i8b) :: n, nelem_sqrt

    nside = -1; nstokes = -1
    if (present(npixtot)) npixtot = -1

    n = get_file_size_dp(file_name)
    if (present(nelem)) nelem = n

    IF (n < 1) RETURN

    nstokes = 1
    IF (DBLE(n/3**4) == DBLE(n)/3.0**4) nstokes = 3

    nelem_sqrt = nint(SQRT(DBLE(n)))

    IF (n /= nelem_sqrt**2) RETURN
    IF (PRESENT(npixtot)) npixtot = nelem_sqrt

    nside = nint(SQRT(DBLE(nelem_sqrt/nstokes/12))) ! npix2nside
    if (nside - nint(SQRT(DBLE(nelem_sqrt/nstokes/12._dp))) /= 0) nside = -1

  END SUBROUTINE get_matrix_size_dp



  FUNCTION wtime()
    REAL(dp) :: wtime
    INTEGER :: count, rate
    !call cpu_time(wtime)
    call system_clock(count, rate)
    wtime = dble(count)/rate
  END FUNCTION wtime



!!$  SUBROUTINE reduce_matrix(matrix, total_hits)
!!$    ! This routine removes zero lines and columns from the matrix for
!!$    ! inversion
!!$    REAL(dp), POINTER :: matrix(:,:)
!!$    INTEGER, POINTER  :: total_hits(:) ! hit mask map
!!$    REAL(dp), POINTER :: matrix_temp(:,:)
!!$    INTEGER           :: npix, sigs, pixel, ierr, hits, &
!!$         a, b, row, col, reduced_row, reduced_col
!!$
!!$    npix = SIZE(total_hits, 1)
!!$    sigs = SIZE(matrix, 1) / npix
!!$    hits  = COUNT(total_hits >= 0)
!!$
!!$    ALLOCATE(matrix_temp(0:sigs*npix-1, 0:sigs*npix-1), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       WRITE (*,*) 'ERROR: unable to allocate in reduce_matrix'
!!$       RETURN
!!$    END IF
!!$    matrix_temp = matrix
!!$    matrix = 0.0
!!$
!!$    reduced_col = 0
!!$    DO col = 0, npix-1
!!$       IF (total_hits(col) == 0) CYCLE
!!$       reduced_row = 0
!!$       DO row = 0, npix-1
!!$         IF (total_hits(row) == 0) CYCLE
!!$         DO a = 0,sigs-1
!!$            DO b = 0,sigs-1
!!$               matrix(reduced_row+a*hits, reduced_col+b*hits) = &
!!$                    matrix_temp(row+a*npix, col+b*npix)
!!$            END DO
!!$         END DO
!!$         reduced_row = reduced_row + 1
!!$       END DO
!!$       reduced_col = reduced_col + 1
!!$    END DO
!!$
!!$    DEALLOCATE(matrix_temp)
!!$
!!$  END SUBROUTINE reduce_matrix
!!$
!!$
!!$
!!$  SUBROUTINE expand_matrix(matrix, total_hits)
!!$    ! This routine restores original matrix size after call to
!!$    ! reduce_matrix()
!!$    REAL(dp), POINTER :: matrix(:,:)
!!$    INTEGER, POINTER  :: total_hits(:)
!!$    REAL(dp), POINTER :: matrix_temp(:,:)
!!$    INTEGER           :: npix, sigs, pixel, ierr, hits, &
!!$         a, b, row, col, reduced_row, reduced_col
!!$
!!$    npix = SIZE(total_hits, 1)
!!$    sigs = SIZE(matrix, 1) / npix
!!$    hits  = COUNT(total_hits >= 0)
!!$
!!$    ALLOCATE(matrix_temp(0:sigs*npix-1, 0:sigs*npix-1), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       WRITE (*,*) 'ERROR: unable to allocate in reduce_matrix'
!!$       RETURN
!!$    END IF
!!$    matrix_temp = matrix
!!$    matrix = 0.0
!!$
!!$    col = -1
!!$    DO reduced_col = 0, hits-1
!!$       DO ! find the real column corresponding to the reduced one
!!$          col = col + 1
!!$          IF (total_hits(col) > 0) EXIT
!!$       END DO
!!$
!!$       row = -1
!!$       DO reduced_row = 0, hits-1
!!$          DO ! find the real row corresponding to the reduced one
!!$             row = row + 1
!!$             IF (total_hits(row) > 0) EXIT
!!$          END DO
!!$         DO a = 0,sigs-1
!!$            DO b = 0,sigs-1
!!$               matrix(row+a*npix, col+b*npix) = &
!!$                    matrix_temp(reduced_row+a*hits, reduced_col+b*hits)
!!$            END DO
!!$         END DO
!!$       END DO
!!$    END DO
!!$
!!$    DEALLOCATE(matrix_temp)
!!$
!!$  END SUBROUTINE expand_matrix
!!$
!!$
!!$
!!$  SUBROUTINE convert_to_interlaced(matrix)
!!$    ! This routine takes a block formed polarization map
!!$    ! and returns it in interlaced form where each
!!$    ! pixel is represented by a 3x3 matrix
!!$    REAL(dp), POINTER :: matrix(:,:)
!!$    REAL(dp), POINTER :: matrix_temp(:,:)
!!$    INTEGER           :: npix, ierr, a, b
!!$
!!$    IF (LBOUND(matrix,1) /= 0 .OR. MODULO(UBOUND(matrix,1)+1,3) /= 0) THEN
!!$       WRITE (*,*) 'ERROR in convert_to_interlaced, matrix does not '//&
!!$            'include polarization'
!!$       RETURN
!!$    END IF
!!$
!!$    npix = (UBOUND(matrix,1)-LBOUND(matrix,1)+1)/3
!!$
!!$    ALLOCATE(matrix_temp(0:npix*3-1, 0:npix*3-1), stat=ierr)
!!$
!!$    IF (ierr/=0) THEN
!!$       WRITE (*,*) 'ERROR : no room to convert into interlaced'
!!$       RETURN
!!$    END IF
!!$
!!$    matrix_temp(:,:) = matrix(:,:)
!!$
!!$    DO a = 0,2
!!$       DO b = 0,2
!!$          matrix(a:3*npix-1:3, b:3*npix-1:3) = &
!!$               matrix_temp(a*npix:(a+1)*npix-1,b*npix:(b+1)*npix-1)
!!$       END DO
!!$    END DO
!!$
!!$    DEALLOCATE(matrix_temp)
!!$
!!$  END SUBROUTINE convert_to_interlaced
!!$
!!$
!!$
!!$  SUBROUTINE convert_to_block_form(matrix)
!!$    ! This routine takes an interlaced polarization map
!!$    ! and returns it in block form where correlations are in
!!$    ! II IQ IU
!!$    ! QU QQ QU
!!$    ! UI UQ UU form
!!$    REAL(dp), POINTER :: matrix(:,:)
!!$    REAL(dp), POINTER :: matrix_temp(:,:)
!!$    INTEGER           :: npix, ierr, a, b
!!$
!!$    IF (LBOUND(matrix,1) /= 0 .OR. MODULO(UBOUND(matrix,1)+1,3) /= 0) THEN
!!$       WRITE (*,*) 'ERROR in convert_to_block_form, matrix does not contain'//&
!!$            ' polarization part'
!!$       RETURN
!!$    END IF
!!$
!!$    npix = (UBOUND(matrix,1)-LBOUND(matrix,1)+1)/3
!!$
!!$    ALLOCATE(matrix_temp(0:npix*3-1, 0:npix*3-1), stat=ierr)
!!$
!!$    IF (ierr/=0) THEN
!!$       WRITE (*,*) 'ERROR : no room to convert into block form'
!!$       RETURN
!!$    END IF
!!$
!!$    matrix_temp(:,:) = matrix(:,:)
!!$
!!$    DO a = 0,2
!!$       DO b = 0,2
!!$          matrix(a*npix:(a+1)*npix-1,b*npix:(b+1)*npix-1) = &
!!$          matrix_temp(a:3*npix-1:3, b:3*npix-1:3)
!!$       END DO
!!$    END DO
!!$
!!$    DEALLOCATE(matrix_temp)
!!$
!!$  END SUBROUTINE convert_to_block_form
!!$
!!$
!!$
!!$  SUBROUTINE save_matrix(matrix, rows, cols, filename, unit)
!!$    REAL(dp), POINTER     :: matrix(:,:)
!!$    INTEGER(dp), OPTIONAL :: rows, cols
!!$    CHARACTER(len=*)      :: filename
!!$    INTEGER, OPTIONAL     :: unit
!!$    INTEGER               :: out_unit=covmat_util_unit, i, j, elems
!!$
!!$    IF (PRESENT(unit)) out_unit = unit
!!$
!!$    elems = SIZE(matrix,1)
!!$
!!$    OPEN(unit=out_unit, file=TRIM(filename), status='replace', recl=elems*15)
!!$    DO i = LBOUND(matrix, 1), UBOUND(matrix, 1)
!!$       DO j = LBOUND(matrix, 2), UBOUND(matrix, 2)
!!$          WRITE (out_unit, '(ES15.5)', advance='no') &
!!$               matrix(i, j)
!!$       END DO
!!$       WRITE (out_unit,*)
!!$    END DO
!!$    CLOSE(out_unit)
!!$
!!$  END SUBROUTINE save_matrix
!!$
!!$
!!$
!!$  SUBROUTINE save_matrix_bin(matrix, rows, cols, filename, unit)
!!$    REAL(dp), POINTER     :: matrix(:,:)
!!$    INTEGER(dp), OPTIONAL :: rows, cols
!!$    CHARACTER(len=*)      :: filename
!!$    INTEGER, OPTIONAL     :: unit
!!$    INTEGER               :: out_unit=covmat_util_unit, i, j
!!$
!!$    IF (PRESENT(unit)) out_unit = unit
!!$
!!$    OPEN(unit=out_unit,file=TRIM(filename),status='replace',form='unformatted')
!!$    WRITE(out_unit) matrix
!!$    CLOSE(out_unit)
!!$
!!$  END SUBROUTINE save_matrix_bin



  SUBROUTINE save_vector(vector, filename, unit)
    REAL(dp), POINTER :: vector(:)
    CHARACTER(len=*) :: filename
    INTEGER, OPTIONAL :: unit
    INTEGER :: out_unit=covmat_util_unit

    IF (PRESENT(unit)) out_unit = unit

    OPEN(unit=out_unit, file=TRIM(filename), status='replace')
    WRITE (out_unit, '(ES20.10)') vector
    CLOSE(out_unit)

  END SUBROUTINE save_vector



!!$  SUBROUTINE save_ivector(vector, filename, unit)
!!$    INTEGER, POINTER      :: vector(:)
!!$    CHARACTER(len=*)      :: filename
!!$    INTEGER, OPTIONAL     :: unit
!!$    INTEGER :: out_unit=covmat_util_unit, i
!!$
!!$    IF (PRESENT(unit)) out_unit = unit
!!$
!!$    OPEN(unit=out_unit, file=TRIM(filename), status='replace')
!!$    WRITE (out_unit, '(I12)') vector
!!$    CLOSE(out_unit)
!!$
!!$  END SUBROUTINE save_ivector



  SUBROUTINE eigen_invert_symmetric_matrix(matrix, limit)
    ! computes the eigenvalues of given matrix and inverts
    ! it using the eigenvalues and vectors
    REAL(dp), POINTER :: matrix(:, :)
    REAL, OPTIONAL :: limit
    REAL :: rcond_limit = 1E-30
    INTEGER(i8b) :: workspace_length
    INTEGER :: ierr, N, i
    REAL(dp), POINTER :: workspace(:), eigenvectors(:,:), eigenvectorsT(:,:), &
         eigenvalues(:)

    IF (PRESENT(limit)) rcond_limit = limit

    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    WRITE (*,'(a,i8)') 'eigen_invert : N == ', N
    ALLOCATE(eigenvalues(N), eigenvectors(N,N), &
         eigenvectorsT(N,N), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') &
            ID, ' : Sorry! no room for eigenvalue inverting matrix'
       RETURN
    END IF

    workspace_length = MAX(5*N, 10000)
    ALLOCATE(workspace(workspace_length), stat=ierr)
    IF (ierr /=0 ) THEN
       WRITE  (*,'(i3,a)') &
            ID, ' : Sorry! no room for workspace in eigenvalue inversion'
       RETURN
    END IF

    ! use dsyev to compute the eigenvalues
    !SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    eigenvectors = matrix
    CALL tic(634)
    CALL DSYEV('V', 'U', N, eigenvectors, N, eigenvalues, &
         workspace, workspace_length, ierr)
    CALL toc('compute eigenvalues', 634)

    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a,i10)') ID, ' : PROBLEM with eigen_invert, info ==',ierr
       RETURN
    END IF

    ! save eigenvalues and eigenvectors
    !CALL save_vector(eigenvalues, 'eigenvalues.dat')
    OPEN(unit=98, file='eigenvalues.dat', status='replace', form='formatted')
    WRITE(98, '(1000000(ES30.20,/))') eigenvalues
    CLOSE(unit=98)
    OPEN(unit=99, file='evecs.dat', status='replace', &
         form='unformatted', access='direct', recl=N**2*8)
    WRITE(99, rec=1) eigenvectors
    CLOSE(unit=99)

    ! construct the inverted matrix
    matrix = 0
    eigenvectors(:,1) = 0
    eigenvectorsT = TRANSPOSE(eigenvectors)
    DO i = 0, N ! ignore worst eigenvalue
       IF ( ABS(eigenvalues(i)/eigenvalues(N)) > rcond_limit ) &
            eigenvectors(:,i) = eigenvalues(i)**(-1)*eigenvectors(:,i)
    END DO
    matrix = MATMUL(eigenvectors, eigenvectorsT)

  END SUBROUTINE eigen_invert_symmetric_matrix



  SUBROUTINE eigen_invert_hermitian_matrix(matrix, limit)
    ! computes the eigenvalues of given matrix and inverts
    ! it using the eigenvalues and vectors
    COMPLEX(dpc), POINTER :: matrix(:, :)
    REAL, OPTIONAL :: limit
    REAL :: rcond_limit = 1E-6
    INTEGER(dp) :: workspace_length
    INTEGER :: ierr, N, i
    COMPLEX(dpc), POINTER :: workspace(:)
    COMPLEX(dpc), POINTER :: eigenvectors(:, :), eigenvectorsT(:, :)
    REAL(dp), POINTER :: eigenvalues(:), workspace2(:)

    IF (PRESENT(limit)) rcond_limit = limit

    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    WRITE (*,'(a,i8)') 'eigen_invert : N == ', N
    ALLOCATE(eigenvalues(N), eigenvectors(N,N), &
         eigenvectorsT(N,N), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') &
            ID, ' : Sorry! no room for eigenvalue inverting matrix'
       RETURN
    END IF

    workspace_length = MAX(5*N, 10000)
    ALLOCATE(workspace(workspace_length), workspace2(3*N-2), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE  (*,'(i3,a)') &
            ID, ' : Sorry! no room for workspace in eigenvalue inversion'
       RETURN
    END IF

    ! use zheev to compute the eigenvalues
    !SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
    eigenvectors = matrix
    CALL tic(634)
    CALL ZHEEV('V', 'U', N, eigenvectors, N, eigenvalues, &
         workspace, workspace_length, workspace2, ierr)
    CALL toc('compute eigenvalues', 634)

    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a,i10)') ID, ' : PROBLEM with eigen_invert, info ==',ierr
       RETURN
    END IF

    ! save eigenvalues
    CALL save_vector(eigenvalues, 'eigenvalues.dat')

    ! save eigenvectors
    open(unit=145, file='eigenvectors.dat', status='replace', &
         form='unformatted', access='direct', recl=N*8*2)
    do i = 1, N
       write (unit=145, rec=i) eigenvectors(:,i)
    end do

    ! construct the inverted matrix
    matrix = 0
    eigenvectorsT = CONJG(TRANSPOSE(eigenvectors))
    DO i = 1, N
       IF (ABS(eigenvalues(i)/eigenvalues(N)) > rcond_limit) THEN
          eigenvectors(:,i) = eigenvalues(i)**(-1)*eigenvectors(:,i)
       ELSE
          eigenvectors( :,i) = 0
          eigenvectorsT(i,:) = 0
          WRITE (*,'(a,i6)') ' Skipping eigenmode #', i
       END IF
    END DO
    matrix = MATMUL(eigenvectors, eigenvectorsT)

  END SUBROUTINE eigen_invert_hermitian_matrix




!!$  SUBROUTINE svd_invert_matrix(matrix, limit)
!!$    ! Computes the pseudoinverse of a matrix
!!$    REAL(dp), POINTER :: matrix(:,:)
!!$    REAL, OPTIONAL    :: limit
!!$    REAL              :: rcond_limit = 1E-3
!!$    INTEGER(dp)       :: row, col, workspace_length
!!$    INTEGER           :: ierr, N, good
!!$    REAL(dp), POINTER :: workspace(:), U(:,:), SIGMA(:), VT(:,:)
!!$    INTEGER, POINTER  :: pivotings(:)
!!$
!!$    IF (PRESENT(limit)) rcond_limit = limit
!!$
!!$    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
!!$    ALLOCATE(SIGMA(N), U(N,N), VT(N,N), stat=ierr)
!!$    IF (ierr /= 0) THEN
!!$       WRITE (*,'(i3,a)') &
!!$            ID, ' : Sorry! no room for SVD inverting matrix'
!!$       RETURN
!!$    END IF
!!$
!!$    workspace_length = MAX(10*N, 10000)
!!$    ALLOCATE(workspace(workspace_length), stat=ierr)
!!$    IF (ierr /=0 ) THEN
!!$       WRITE  (*,'(i3,a)') &
!!$            ID, ' : Sorry! no room for workspace in svd inversion'
!!$       RETURN
!!$    END IF
!!$
!!$    ! decompose:
!!$    ! This routine returns the singular value decomposition
!!$    ! of matrix A : A = U * SIGMA * transpose(V), where
!!$    ! U and V are orthogonal matrices and SIGMA is a diagonal matrix
!!$    ! with the singular values in decreasing order on the diagonal
!!$    !
!!$    !SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
!!$    !                   WORK, LWORK, INFO )
!!$    CALL DGESVD('A', 'A', N, N, matrix, N, SIGMA, U, N, VT, N, &
!!$         workspace, workspace_length, ierr)
!!$
!!$    IF (ierr /= 0) THEN
!!$       WRITE (*,'(i3,a)') ID, ' : PROBLEM with SVD'
!!$       RETURN
!!$    END IF
!!$
!!$    CALL save_vector(SIGMA, 'singular_values.dat')
!!$
!!$    ! Check how many good singular values we have. Omit bad ones from
!!$    ! the pseudoinverse
!!$    good = 2
!!$    DO
!!$       IF (good > N) EXIT
!!$       IF (SIGMA(good)/SIGMA(1) < rcond_limit) EXIT
!!$       good = good + 1
!!$    END DO
!!$    good = good - 1
!!$
!!$    !good = 189 ! remove after testing
!!$
!!$    WRITE (*,'(i3,a,i6,a,i6)') ID, ' : Number of good singular values : ', &
!!$         good, '/', N
!!$
!!$    DO col = 1,good
!!$       U(:, col) = U(:, col)/SIGMA(col)
!!$    END DO
!!$
!!$    ! The pseudoinverse is now:
!!$    ! A^-1 = V * SIGMA^inv * transpose(U), where
!!$    ! SIGMA^inv is the original diagonal SIGMA-matrix, with the nonzero
!!$    ! singular values inverted
!!$
!!$    matrix = TRANSPOSE(MATMUL(U(:,1:good),VT(1:good,:)))
!!$
!!$    DEALLOCATE(SIGMA, U, VT, workspace)
!!$
!!$  END SUBROUTINE svd_invert_matrix
!!$
!!$
!!$
!!$  SUBROUTINE invert_matrix(matrix)
!!$    REAL(dp), POINTER :: matrix(:,:)
!!$    INTEGER(dp)       :: row, col, workspace_length, N
!!$    INTEGER           :: ierr
!!$    REAL(dp), POINTER :: workspace(:)
!!$    INTEGER, POINTER  :: pivotings(:)
!!$
!!$    ! factorize, dsytrf assumes covmat_inv to be in upper triangular form
!!$    ! SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!!$    ! LDA=leading dimension, tells the routine how many elements there
!!$    ! are between two consecutive columns on the same row
!!$    ierr = 0
!!$    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
!!$    workspace_length = N**2 * 5
!!$    ALLOCATE(workspace(workspace_length), pivotings(N), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       WRITE (*,'(i3,a)') &
!!$            ID, ' : Sorry, no room for inverting matrix'
!!$    END IF
!!$
!!$    ! Factorize
!!$    CALL DSYTRF('upper', N, matrix, N, pivotings, workspace, &
!!$         workspace_length, ierr)
!!$    IF (ierr /= 0) THEN
!!$       WRITE (*,'(i3,a)') &
!!$            ID,' : PROBLEM with factorizing the covariance matrix'
!!$    END IF
!!$
!!$    ! Invert
!!$    CALL DSYTRI('upper', N, matrix, N, pivotings, workspace, &
!!$         workspace_length, ierr)
!!$    IF (ierr /= 0) THEN
!!$       WRITE (*,'(i3,a)') &
!!$            ID, ' : PROBLEM with inverting the covariance matrix'
!!$    END IF
!!$
!!$    ! Finally just copy the upper triangle to lower half as well
!!$    DO col = 1, N-1
!!$       DO row = 0, col-1
!!$          matrix(col, row) = matrix(row, col)
!!$       END DO
!!$    END DO
!!$
!!$    DEALLOCATE(workspace, pivotings)
!!$
!!$  END SUBROUTINE invert_matrix



  SUBROUTINE tic(index)
    INTEGER, OPTIONAL :: index

    IF (PRESENT(index)) THEN
       tic_times(index) = wtime()
    ELSE
       time_start = wtime()
    END IF

  END SUBROUTINE tic



  SUBROUTINE toc(label, index)
    CHARACTER(len=*), OPTIONAL :: label
    INTEGER, OPTIONAL          :: index
    CHARACTER(len=512)         :: middle_string
    REAL(dp)                   :: elapsed, time_now

    time_now = wtime()

    IF (PRESENT(index)) THEN
       elapsed = time_now - tic_times(index)
    ELSE
       elapsed = time_now - time_start
    END IF

    IF (PRESENT(label)) THEN
       middle_string = ' : '//TRIM(label)//' completed in '
    ELSE
       middle_string = ' : elapsed time : '
    END IF

    WRITE (*,'(i3,a,f8.2,a)') &
         ID, TRIM(middle_string), elapsed, ' s'

  END SUBROUTINE toc


END MODULE covmat_util
