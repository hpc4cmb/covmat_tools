! Collection of routines to use the scalapack routines
!
! 2010-04-27 Reijo Keskitalo


MODULE scalapack_tools

  USE healpix_types, only : sp, dp, i4b, i8b, filenamelen
  USE covmat_util, only : get_file_size_dp, get_matrix_size_dp

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE scalapack_env
     INTEGER(i4b) :: id, context, nprocs, nprow, npcol, myrow, mycol, BLOCK
  END TYPE scalapack_env

  TYPE scalapack_matrix_dp
     INTEGER(i4b) :: nrow_tot, ncol_tot, nrow, ncol, ld, desc(9), filetype
     INTEGER(i4b) :: ndata
     REAL(dp), POINTER :: DATA(:)=>null()
  END TYPE scalapack_matrix_dp

  EXTERNAL numroc
  INTEGER(i4b) :: numroc

  private :: sp, dp, i4b, i8b, filenamelen


CONTAINS



  SUBROUTINE init_scalapack(spenv, blocksize, nprow_in, npcol_in)

    IMPLICIT NONE

    TYPE(scalapack_env) :: spenv
    INTEGER(i4b), OPTIONAL :: blocksize, nprow_in, npcol_in

    INTEGER(i4b) :: nprow, npcol

    CALL BLACS_PINFO(spenv%id, spenv%nprocs)

    if (present(nprow_in)) then
       ! user has specified the number of rows for the processor matrix
       nprow = nprow_in
       npcol = spenv%nprocs/nprow
       if (nprow*npcol /= spenv%nprocs .and. spenv%id == 0) &
            write (*,'(" WARNING: user-specified processor matrix is sub-optimal")')
    else if (present(npcol_in)) then
       ! user has specified the number of columns for the processor matrix
       npcol = npcol_in
       nprow = spenv%nprocs/npcol
       if (nprow*npcol /= spenv%nprocs .and. spenv%id == 0) &
            write (*,'(" WARNING: user-specified processor matrix is sub-optimal")')
    else
       ! choose the processor grid to be maximally symmetric
       DO nprow = INT(SQRT(DBLE(spenv%nprocs))),1,-1
          npcol = spenv%nprocs/nprow
          IF (npcol * nprow == spenv%nprocs) EXIT
       END DO
    end if
    spenv%nprow = nprow
    spenv%npcol = npcol

    spenv%block = 32
    IF (PRESENT(blocksize)) spenv%block = blocksize

    IF (spenv%id == 0) &
         WRITE(*,'(/," Organizing processes into a ",i0,"x",i0,&
         & " grid, blocksize == ",i0)') &
         spenv%nprow, spenv%npcol, spenv%block

    ! get default system context
    CALL BLACS_GET(0, 0, spenv%context)

    ! form a grid of tasks
    CALL BLACS_GRIDINIT(spenv%context, 'column-major', &
         spenv%nprow, spenv%npcol)

    ! compute my location on the grid
    CALL BLACS_GRIDINFO(spenv%context, spenv%nprow, spenv%npcol, &
         spenv%myrow, spenv%mycol)

#ifdef INTEL_MKL
    CALL mkl_set_num_threads(1)
#endif

  END SUBROUTINE init_scalapack



  SUBROUTINE finalize_scalapack()

    IMPLICIT NONE

    ! frees memory held by all contexts
    ! but appears to crash with the default system context
    ! CALL BLACS_EXIT(1)

  end SUBROUTINE finalize_scalapack



  SUBROUTINE free_scalapack_matrix_dp(matrix)

    IMPLICIT NONE

    TYPE(scalapack_matrix_dp) :: matrix

    if (associated(matrix%data)) deallocate(matrix%data)

  END SUBROUTINE free_scalapack_matrix_dp



  SUBROUTINE initialize_scalapack_matrix_dp(spenv, matrix, nrow_tot, ncol_tot)

    IMPLICIT NONE

    TYPE(scalapack_env) :: spenv
    TYPE(scalapack_matrix_dp) :: matrix
    INTEGER(i4b), OPTIONAL :: nrow_tot, ncol_tot

    INTEGER(i4b) :: ierr, info
    INTEGER(i4b) :: ndims, false_id
    INTEGER(i4b), POINTER :: gsizes(:), distribs(:), dargs(:), psizes(:)

    ! if user did not supply dimensions, they are assumed to be stored
    ! in the matrix structure
    IF (PRESENT(nrow_tot)) matrix%nrow_tot = nrow_tot
    IF (PRESENT(ncol_tot)) matrix%ncol_tot = ncol_tot

    ! consistency checks
    IF (matrix%nrow_tot <= 0 .or. matrix%ncol_tot <= 0) THEN
       WRITE (*,*) 'ERROR in initialize_scalapack_matrix: ' // &
            'matrix dimensions must be positive'
       RETURN
    END IF


    matrix%nrow = NUMROC(matrix%nrow_tot, spenv%block, spenv%myrow, 0, &
         spenv%nprow)
    matrix%ncol = NUMROC(matrix%ncol_tot, spenv%block, spenv%mycol, 0, &
         spenv%npcol)

    matrix%ndata = matrix%nrow * matrix%ncol

    ALLOCATE(matrix%data(matrix%ndata), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,*) 'ERROR in initialize_scalapack_matrix: ' // &
            'no room for matrix!'
       RETURN
    END IF

    matrix%data = 0.0_dp


    ! create the descriptor
    CALL DESCINIT(matrix%desc, matrix%nrow_tot, matrix%ncol_tot, &
         spenv%block, spenv%block, 0, 0, spenv%context, matrix%nrow, info)

    ! create an mpi_datatype for the file
    ndims = 2
    ALLOCATE(gsizes(ndims), distribs(ndims), dargs(ndims), psizes(ndims))
    gsizes = (/ matrix%nrow_tot, matrix%ncol_tot /) ! global size of the matrix
    distribs = MPI_DISTRIBUTE_CYCLIC ! distribution for each dimension
    dargs = spenv%block ! additional distribution parameters
    psizes = (/ spenv%nprow, spenv%npcol /) ! size of the process grid

    ! darray orders the processes row-major whereas scalapack assumes
    ! them to be column-major. We must transpose the processor grid
    false_id = spenv%myrow*spenv%npcol + spenv%mycol
    CALL MPI_TYPE_CREATE_DARRAY(spenv%nprocs, false_id, ndims, &
         gsizes, distribs, dargs, psizes, MPI_ORDER_FORTRAN, &
         MPI_REAL8, matrix%filetype, ierr)
    IF (ierr /= MPI_SUCCESS) STOP 'ERROR in mpi_type_create_darray'

    CALL MPI_TYPE_COMMIT(matrix%filetype, ierr)
    IF (ierr /= MPI_SUCCESS) STOP 'ERROR in mpi_type_commit'

    CALL MPI_BARRIER(mpi_comm_world, ierr)

  END SUBROUTINE initialize_scalapack_matrix_dp



  SUBROUTINE read_scalapack_matrix_dp(spenv, matrix, file_matrix, &
       nrow_tot, ncol_tot)

    IMPLICIT NONE

    TYPE(scalapack_env) :: spenv
    TYPE(scalapack_matrix_dp) :: matrix
    CHARACTER(len=*) :: file_matrix
    INTEGER(i8b), OPTIONAL :: nrow_tot, ncol_tot

    LOGICAL :: there
    INTEGER(i4b) :: ierr
    INTEGER(i4b) :: status(MPI_STATUS_SIZE), filemode, fileinfo, handle_matrix
    INTEGER(MPI_OFFSET_KIND) :: fileoffset

    IF (spenv%id == 0) WRITE (*,'(/," Reading ",a,/)') TRIM(file_matrix)

    ! if user did not supply dimensions, they are assumed to be stored
    ! in the matrix structure
    IF (PRESENT(nrow_tot)) matrix%nrow_tot = int(nrow_tot, i4b)
    IF (PRESENT(ncol_tot)) matrix%ncol_tot = int(ncol_tot, i4b)
    
    IF (spenv%id == 0) then
       WRITE (*,'(" nrow = ",i0,", ncol = ", i0,/)') &
            matrix%nrow_tot, matrix%ncol_tot
    end IF

    CALL initialize_scalapack_matrix_dp(spenv, matrix)

    ! consistency checks
    INQUIRE(file=TRIM(file_matrix), exist=there)
    IF (.NOT. there) THEN
       IF (spenv%id == 0) &
            WRITE (*,*) 'ERROR in read_scalapack_matrix_dp: ' // &
            'matrix does not exist!'
       RETURN
    END IF

    ! open the mpi-io file view
    CALL mpi_info_create(fileinfo, ierr)
    filemode = MPI_MODE_RDONLY
    CALL mpi_file_open(mpi_comm_world, TRIM(file_matrix), filemode, fileinfo, &
         handle_matrix, ierr)
    IF (ierr /= 0) THEN
       IF (spenv%id == 0) &
            WRITE (*,*) 'ERROR in read_scalapack_matrix_dp: ' // &
            'Unable to open file'
       RETURN
    END IF

    fileoffset = 0
    CALL mpi_file_set_view(handle_matrix, fileoffset, MPI_REAL8, &
         matrix%filetype, 'native', fileinfo, ierr)
    IF (ierr /= 0) THEN
       IF (spenv%id == 0) &
            WRITE (*,*) 'ERROR in read_scalapack_matrix_dp: ' // &
            'Unable to establish view to matrix file'
       RETURN
    END IF


    CALL mpi_file_read_all(handle_matrix, matrix%data, matrix%ndata, &
         MPI_REAL8, status, ierr)
    IF (ierr /= MPI_SUCCESS) STOP 'ERROR in mpi_file_read_all'


    ! finalize
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_file_close(handle_matrix, ierr)

  END SUBROUTINE read_scalapack_matrix_dp



  SUBROUTINE write_scalapack_matrix_dp(spenv, matrix, file_matrix)

    IMPLICIT NONE

    TYPE(scalapack_env) :: spenv
    TYPE(scalapack_matrix_dp) :: matrix
    CHARACTER(len=*) :: file_matrix

    INTEGER(i4b) :: ierr
    INTEGER(i4b) :: status(MPI_STATUS_SIZE), filemode, fileinfo, handle_matrix
    INTEGER(MPI_OFFSET_KIND) :: fileoffset

    IF (spenv%id == 0) WRITE (*,'(/," Writing ",a,/)') TRIM(file_matrix)

    ! open the mpi-io file view
    CALL mpi_info_create(fileinfo, ierr)
    filemode = IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE)
    CALL mpi_file_open(mpi_comm_world, file_matrix, filemode, fileinfo, &
         handle_matrix, ierr)
    IF (ierr /= 0) THEN
       IF (spenv%id == 0) &
            WRITE (*,*) 'ERROR in write_scalapack_matrix_dp: ' // &
            'Unable to open file'
       RETURN
    END IF

    fileoffset = 0
    CALL mpi_file_set_view(handle_matrix, fileoffset, MPI_REAL8, &
         matrix%filetype, 'native', fileinfo, ierr)
    IF (ierr /= 0) THEN
       IF (spenv%id == 0) &
            WRITE (*,*) 'ERROR in write_scalapack_matrix_dp: ' // &
            'Unable to establish view to matrix file'
       RETURN
    END IF

    ! write the file
    CALL mpi_file_write_all(handle_matrix, matrix%data, matrix%ndata, &
         MPI_REAL8, status, ierr)
    IF (ierr /= MPI_SUCCESS) STOP 'ERROR in mpi_file_write_all'

    ! finalize
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_file_close(handle_matrix, ierr)

  END SUBROUTINE write_scalapack_matrix_dp



  SUBROUTINE eigendecomp_scalapack_matrix_dp(spenv, matrix, evals, evecs)
    
    IMPLICIT NONE
    
    TYPE(scalapack_env) :: spenv
    TYPE(scalapack_matrix_dp) :: matrix, evecs
    REAL(dp), POINTER :: evals(:)
    
    REAL(dp), POINTER :: workspace_dp(:)
    INTEGER(i4b), POINTER :: workspace_i4b(:)
    INTEGER(i4b) :: info, ierr
    INTEGER(i8b) :: trilwmin, len_workspace_dp, len_workspace_i4b

    if (spenv%id == 0) write (*,'(/," Decomposing matrix",/)')

    ! consistency check
    IF (matrix%nrow_tot /= matrix%ncol_tot .OR. &
         matrix%nrow_tot <= 0 .or.  matrix%ncol_tot <= 0) THEN
       WRITE (*,'(i4," : ",a,i0,a,i0)') spenv%id, &
            'ERROR in eigendecomp_scalapack_matrix_dp: matrix dimensions ' // &
            'are not symmetric: ', matrix%nrow_tot, ' x ', matrix%ncol_tot
       RETURN
    END IF

    ! initialize the output matrix
    CALL initialize_scalapack_matrix_dp(spenv, evecs, matrix%nrow_tot, &
         matrix%ncol_tot)

    ALLOCATE(evals(matrix%nrow_tot), stat=ierr)
    IF (ierr /= 0) THEN
       IF (spenv%id == 0) &
            WRITE (*,*) 'ERROR in eigendecomp_scalapack_matrix_dp: ' // &
            'No room for eigenvalues'
       RETURN
    END IF

    ! see http://www.netlib.org/scalapack/double/pdsyevd.f
    !  SUBROUTINE PDSYEVD( JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
    ! $                    DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )


    ! query workspace size
    !          LWORK >= MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
    !          TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
    !          NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
    !          NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
    len_workspace_dp = 1 + 6*matrix%nrow_tot + 2*matrix%ncol*matrix%nrow
    trilwmin = 3*matrix%nrow_tot + &
         nint(MAXVAL(DBLE((/ spenv%block*(matrix%nrow+1), 3*spenv%block /))))
    IF (trilwmin > len_workspace_dp) len_workspace_dp = trilwmin
    len_workspace_dp = len_workspace_dp + 2*matrix%nrow_tot
    IF (len_workspace_dp < 200._dp*2**20/8) len_workspace_dp = 200*2**20/8

    len_workspace_i4b = 7*matrix%nrow_tot + 8*spenv%npcol + 2
    IF (len_workspace_i4b < 50._dp*2**20/4) len_workspace_i4b = 50*2**20/4

    ALLOCATE(workspace_dp(len_workspace_dp), workspace_i4b(len_workspace_i4b), &
         stat=ierr)
    IF (ierr /= 0) THEN
       IF (spenv%id == 0) &
            WRITE (*,*) 'ERROR in eigendecomp_scalapack_matrix_dp: ' // &
            'No room for workspace'
       RETURN
    END IF

    ! perform the eigenvalue decomposition
    CALL pdsyevd('V', 'U', matrix%nrow_tot, matrix%data, 1, 1, matrix%desc, &
         evals, evecs%data, 1, 1, evecs%desc, &
         workspace_dp, len_workspace_dp, &
         workspace_i4b, len_workspace_i4b, info)

    ! check the error flag
    IF (info /= 0) THEN
       IF (spenv%id == 0) &
            WRITE (*,'(a,i0)') 'ERROR in eigendecomp_scalapack_matrix_dp: ' // &
            'pdsyevd failed with info == ', info
       RETURN
    END IF

    CALL mpi_barrier(mpi_comm_world, ierr)

  END SUBROUTINE eigendecomp_scalapack_matrix_dp



  SUBROUTINE multiply_scalapack_matrix_dp(spenv, matrix_in1, matrix_in2, &
       matrix_out, ktranspose)
    
    IMPLICIT NONE
    
    TYPE(scalapack_env) :: spenv
    TYPE(scalapack_matrix_dp) :: matrix_in1, matrix_in2, matrix_out
    LOGICAL, OPTIONAL :: ktranspose(2)
    CHARACTER(len=80) :: string_transpose(2)
    INTEGER(i4b) :: nrow_out, ncol_out, ncol_in

    if (spenv%id == 0) write (*,'(/," Multiplying matrices",/)')

    string_transpose = 'No transpose'
    IF (PRESENT(ktranspose)) WHERE(ktranspose) string_transpose='Transpose'
    
    ! get the dimensions of the output matrix
    IF (INDEX(string_transpose(1), 'No') > 0) THEN
       nrow_out = matrix_in1%nrow_tot
       ncol_in = matrix_in1%ncol_tot
    ELSE
       nrow_out = matrix_in1%ncol_tot
       ncol_in = matrix_in1%nrow_tot
    END IF

    IF (INDEX(string_transpose(2), 'No') > 0) THEN
       ncol_out = matrix_in2%ncol_tot
       IF (ncol_in /= matrix_in2%nrow_tot) THEN
          WRITE (*,*) 'ERROR in multiply_scalapack_matrix_dp: ' // &
               'inner matrix dimensions do not match'
          RETURN
       END IF
    ELSE
       ncol_out = matrix_in2%nrow_tot
       IF (ncol_in /= matrix_in2%ncol_tot) THEN
          WRITE (*,*) 'ERROR in multiply_scalapack_matrix_dp: ' // &
               'inner matrix dimensions do not match'
          RETURN
       END IF
    END IF

    ! initialize the output matrix
    CALL initialize_scalapack_matrix_dp(spenv, matrix_out, nrow_out, ncol_out)

    ! Multiply
    CALL pdgemm(string_transpose(1), string_transpose(2), &
         nrow_out, ncol_out, ncol_in, &
         1.0_dp, matrix_in1%data, 1, 1, matrix_in1%desc, &
         matrix_in2%data, 1, 1, matrix_in2%desc, &
         0.0_dp, matrix_out%data, 1, 1, matrix_out%desc)

  END SUBROUTINE multiply_scalapack_matrix_dp



  SUBROUTINE symmetrize_scalapack_matrix_dp(spenv, matrix_in)
    
    IMPLICIT NONE
    
    TYPE(scalapack_env) :: spenv
    TYPE(scalapack_matrix_dp) :: matrix_in

    INTEGER(i4b) :: nrow, ncol, npix

    if (spenv%id == 0) write (*,'(/," Symmetrizing matrix",/)')

    ! consistency check
    nrow = matrix_in%nrow_tot
    ncol = matrix_in%ncol_tot
    if (nrow /= ncol .or. nrow < 1) then
       write (*,'("ERROR: cannot symmetrize non-square matrix")')
       return
    end if
    npix = nrow

    call pdtran(npix, npix, &
         0.5_dp, matrix_in%data, 1, 1, matrix_in%desc, &
         0.5_dp, matrix_in%data, 1, 1, matrix_in%desc)

  END SUBROUTINE symmetrize_scalapack_matrix_dp
  


END MODULE scalapack_tools
