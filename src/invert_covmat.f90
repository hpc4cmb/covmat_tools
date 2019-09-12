! This code reads madam pixel-pixel inverse covariance matrix
! and inverts it.
!
! Revision history:
!   2012-10-25 : nside and nstokes are now measured from the file size,
!                eigenlimit is optional

PROGRAM invert_madam_covmat

  USE healpix_types
  USE extension, only   : nArguments, getArgument
  USE covmat_util, only : tic, toc, get_matrix_size_dp

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfile
  INTEGER(i4b) :: nstokes, nside
  integer(i8b) :: npix
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:,:)
  LOGICAL :: there
  INTEGER(i4b) :: ierr, rec_len
  real(dp) :: limit

  IF (nArguments() /= 1 .and. nArguments() /= 2)  THEN
     STOP 'Usage: invert_covmat <covmat> [<limit>]'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     if (nArguments() == 2) then
        CALL getArgument(2, argument)
        read(argument, *, iostat=ierr) limit
        if ( ierr /= 0 ) then
           print *,'Failed to parse limit from ' // trim(argument)
           stop
        end if
        write (*,*) 'Using threshold eval > evalmax * ',limit
     end if
  END IF

  call get_matrix_size_dp(covmat_file, nside, nstokes, npix)
  WRITE (*,'(a,i8)') '  Nside == ', nside
  WRITE (*,'(a,i8)') '   Npix == ', npix
  WRITE (*,'(a,i8)') 'Nstokes == ', nstokes
  WRITE (*,'(a,i8)') 'Npixtot == ', npix*nstokes

  ALLOCATE(covmat(0:npix*nstokes-1, 0:npix*nstokes-1), stat=ierr)
  IF (ierr/=0) STOP 'no room for covmat'

  call tic
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old',&
       form='unformatted', access='direct', recl=rec_len)
  READ (covmat_unit,rec=1) covmat
  CLOSE(unit=covmat_unit)
  call toc('Read matrix')

  call tic
  call eigen_invert_symmetric_matrix(covmat, covmat_file, limit)
  call toc('Invert matrix')

  call tic
  outfile = TRIM(ADJUSTL(TRIM(covmat_file)//'.inverse'))
  write (*,*) ' Writing to: ', TRIM(outfile)
  OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(outfile)), status='replace', &
       form='unformatted', access='direct', recl=rec_len)
  WRITE (covmat_unit, rec=1) covmat
  CLOSE(unit=covmat_unit)
  call toc('Write matrix')

  deallocate(covmat)

CONTAINS

  SUBROUTINE eigen_invert_symmetric_matrix(matrix, covmat_file, limit)
    ! computes the eigenvalues of given matrix and inverts
    ! it using the eigenvalues and vectors
    REAL(dp), POINTER :: matrix(:,:)
    CHARACTER(len=*)  :: covmat_file
    REAL(dp), OPTIONAL :: limit
    REAL(dp) :: rcond_limit = 1E-30
    INTEGER(i8b) :: workspace_length
    INTEGER(i4b) :: ierr, N, i
    REAL(dp), POINTER :: workspace(:), eigenvectors(:,:), eigenvectorsT(:,:), &
         eigenvalues(:)

    IF (PRESENT(limit)) rcond_limit = limit

    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    WRITE (*,'(a,i8)') 'eigen_invert : N == ', N
    ALLOCATE(eigenvalues(N), eigenvectors(N,N), &
         eigenvectorsT(N,N), stat=ierr)
    IF (ierr /= 0) THEN
       stop ' Sorry! no room for eigenvalue inverting matrix'
    END IF

    workspace_length = MAX(5*N, 10000)
    ALLOCATE(workspace(workspace_length), stat=ierr)
    IF (ierr /=0 ) THEN
       stop ' Sorry! no room for workspace in eigenvalue inversion'
    END IF

    ! use dsyev to compute the eigenvalues
    !SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    eigenvectors = matrix
    CALL tic(634)
    CALL DSYEV('V', 'U', N, eigenvectors, N, eigenvalues, &
         workspace, workspace_length, ierr)
    CALL toc('compute eigenvalues', 634)

    IF (ierr /= 0) THEN
       WRITE (*,'(a,i10)') ' PROBLEM with eigen_invert, info ==',ierr
       stop
    END IF

    ! save eigenvalues and eigenvectors
    OPEN(unit=98, file=trim(covmat_file)//'.eigenvalues', status='replace', &
         form='formatted')
    WRITE(98, '(1000000(ES30.20,/))') eigenvalues
    CLOSE(unit=98)
    OPEN(unit=99, file=trim(covmat_file)//'-evecs', status='replace', &
         form='unformatted', access='direct', recl=rec_len)
    WRITE(99, rec=1) eigenvectors
    CLOSE(unit=99)

    ! invert the eigenvalues for recomposing
    where(eigenvalues < maxval(eigenvalues)*rcond_limit) eigenvalues = 0
    where(eigenvalues /= 0) eigenvalues = 1.0 / sqrt(eigenvalues)
    OPEN(unit=98, file=trim(covmat_file)//'.eigen', status='replace', &
         form='formatted')
    WRITE(98, '(1000000(ES30.20,/))') eigenvalues
    CLOSE(unit=98)

    ! construct the inverted matrix    
    matrix = 0
    eigenvectorsT = TRANSPOSE(eigenvectors)
    DO i = 1, N
       eigenvectors(:,i) = eigenvalues(i)*eigenvectors(:,i)
    END DO
    matrix = MATMUL(eigenvectors, eigenvectorsT)

    deallocate(eigenvalues, eigenvectors, eigenvectorsT, workspace)

  END SUBROUTINE eigen_invert_symmetric_matrix


END PROGRAM invert_madam_covmat
