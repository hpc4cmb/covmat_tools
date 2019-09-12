! This code reads madam pixel-pixel inverse covariance matrix
! and inverts it.

PROGRAM invert_madam_covmat

  USE healpix_types
  USE extension, only   : nArguments, getArgument
  USE covmat_util, only : tic, toc

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfile
  INTEGER(I8B) :: npix, nstokes
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmat(:, :)
  LOGICAL :: there, interlaced=.false.
  INTEGER(i4b) :: ierr, rec_len

  IF (nArguments() < 3)  THEN
     STOP 'Usage: invert_covmat <covmat> <npix> <nstokes> [<--interlaced>]'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) npix
     WRITE (*,'(a,i4)') 'Npix   == ', npix

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     if (nArguments() == 4) then
        CALL getArgument(4, argument)
        if (index(argument, '-i') > 0) then
           WRITE (*,'(a)') 'Assuming matrix is interlaced'
           interlaced = .true.
        else
           WRITE (*,'(a)') 'Unrecognized option : '//trim(argument)
        end if
     end if
  END IF

  ALLOCATE(covmat(npix*nstokes, npix*nstokes), stat=ierr)
  IF (ierr/=0) STOP 'no room for covmat'

  call tic
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old',&
       form='unformatted', access='direct', recl=rec_len)
  READ (covmat_unit,rec=1) covmat
  CLOSE(unit=covmat_unit)
  call toc('Read matrix')

  covmat = 0.5*(covmat+transpose(covmat))

  call tic
  call eigen_invert_symmetric_matrix(covmat, 1E-6)
  call toc('Invert matrix')

  ! remove the Stokes I offset
  call tic
  if (interlaced) then
     covmat(1:nstokes*npix:nstokes,1:nstokes*npix:nstokes) = &
          covmat(1:nstokes*npix:nstokes,1:nstokes*npix:nstokes) &
          - sum(covmat(1:nstokes*npix:nstokes,1:nstokes*npix:nstokes))/npix**2
  else
     covmat(1:npix,1:npix) = &
          covmat(1:npix,1:npix) - sum(covmat(1:npix,1:npix))/npix**2
  end if
  call toc('Remove I offset')  

  call tic
  outfile = TRIM(ADJUSTL(TRIM(covmat_file)//'.inverse'))
  write (*,*) ' Writing to: ', TRIM(outfile)
  OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(outfile)), status='replace', &
       form='unformatted', access='direct', recl=rec_len)
  WRITE (covmat_unit, rec=1) covmat
  CLOSE(unit=covmat_unit)
  call toc('Write matrix')

  !deallocate(covmat)


CONTAINS


  SUBROUTINE eigen_invert_symmetric_matrix(matrix, limit)
    ! computes the eigenvalues of given matrix and inverts
    ! it using the eigenvalues and vectors
    REAL(dp), POINTER :: matrix(:,:)
    REAL, OPTIONAL    :: limit
    REAL :: rcond_limit = 1E-30
    INTEGER(dp) :: workspace_length
    INTEGER :: ierr, N, i
    REAL(dp), POINTER :: workspace(:), eigenvectors(:,:), eigenvectorsT(:,:),&
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
    OPEN(unit=98, file=trim(covmat_file)//'.eigen', status='replace', &
         form='formatted')
    WRITE(98, '(1000000(ES30.20,/))') eigenvalues
    CLOSE(unit=98)
    OPEN(unit=99, file=trim(covmat_file)//'-evecs', status='replace', &
         form='unformatted', access='direct', recl=rec_len)
    WRITE(99, rec=1) eigenvectors
    CLOSE(unit=99)

    ! construct the inverted matrix    
    matrix = 0
    !eigenvectors(:,1) = 0.0 ! eliminate I monopole, the worst eigenvector
    eigenvectorsT = TRANSPOSE(eigenvectors)
    DO i = 1, N ! ignore worst eigenvalue
       IF ( eigenvalues(i)/eigenvalues(N) > rcond_limit ) &
            eigenvectors(:,i) = eigenvalues(i)**(-1)*eigenvectors(:,i)
    END DO
    matrix = MATMUL(eigenvectors, eigenvectorsT)

  END SUBROUTINE eigen_invert_symmetric_matrix


END PROGRAM invert_madam_covmat
