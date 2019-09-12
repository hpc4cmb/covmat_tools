!> @file
!! Simple interface to PBLAS routine PDGEMM to
!! multiply two block cyclic distributed double precision matrices
!
!! @date 2008-10-26
!! @author Reijo Keskitalo

!> Usage: covmatmul <covmat1> <covmat2> <blocksize> [<out>]
!!
!! @param covmat1 First covmat to multiply
!! @param covmat2 Second covmat to multiply
!! @param blocksize Scalapack blocking factor, typically 32
!! @param out Optional name of the output matrix, default is out.dat
PROGRAM covmatmul

  USE healpix_types
  USE extension, ONLY : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc
  USE scalapack_tools
  USE pix_tools, ONLY : nside2npix

  IMPLICIT NONE

  TYPE(scalapack_env) :: spenv
  TYPE(scalapack_matrix_dp) :: covmat1, covmat2, covmat3
  CHARACTER(len=filenamelen) :: argument, covmat_file, covmat_file2, outfile
  INTEGER, PARAMETER :: covmat_unit=55
  LOGICAL :: there
  INTEGER(i4b) :: ierr, id, ntasks, blocksize
  INTEGER(i8b) :: npix, npix2

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (nArguments() < 3)  THEN
     IF (id == 0) WRITE(*,'(/,a,/)') &
          '  Usage: covmatmul '// &
          '<covmat1> <covmat2> <blocksize> [<out>]'
     CALL mpi_finalize(ierr)
     STOP
  ELSE IF (id == 0) THEN
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     covmat_file2 = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file2), exist=there)
     IF (.NOT. there) STOP 'covmat_file2 does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file2)

     CALL getArgument(3, argument)
     READ (argument, *, iostat=ierr) blocksize
     IF (ierr /= 0) STOP 'Unable to parse block size'
     WRITE (*,'(a,i4)') 'block   == ', blocksize

     outfile = 'out.dat'
     IF (nArguments() > 3) THEN
        CALL getArgument(4, argument)
        outfile = TRIM(argument)
     END IF
     WRITE (*,*) ' writing to ' // TRIM(outfile)

  END IF

  CALL mpi_bcast(covmat_file, filenamelen,mpi_character,0,mpi_comm_world,ierr)
  CALL mpi_bcast(covmat_file2,filenamelen,mpi_character,0,mpi_comm_world,ierr)
  CALL mpi_bcast(blocksize,   1,          mpi_integer,  0,mpi_comm_world,ierr)
  CALL mpi_bcast(outfile,     filenamelen,mpi_character,0,mpi_comm_world,ierr)

  npix = NINT(SQRT(DBLE(get_file_size_dp(covmat_file))))
  npix2 = NINT(SQRT(DBLE(get_file_size_dp(covmat_file2))))
  IF (id == 0) WRITE (*,'(" Npix  == ",i4)') npix
  IF (npix /= npix2) STOP 'Matrix sizes do not match'

  CALL init_scalapack(spenv, blocksize)

  ! read in the matrices
  IF (ID == 0) CALL tic
  CALL read_scalapack_matrix_dp(spenv, covmat1, covmat_file, npix, npix)
  CALL read_scalapack_matrix_dp(spenv, covmat2, covmat_file2, npix, npix) 
     
  CALL mpi_barrier(mpi_comm_world, ierr)

  IF (ID == 0) CALL toc('Read matrices')

  ! Multiply
  IF (ID == 0) CALL tic
  CALL multiply_scalapack_matrix_dp(spenv, covmat1, covmat2, covmat3)
  IF (ID == 0) CALL toc('Multiplication')
     
  ! save the results
  IF (ID == 0) CALL tic  
  CALL write_scalapack_matrix_dp(spenv, covmat3, outfile)
  CALL mpi_barrier(mpi_comm_world, ierr)
  IF (ID == 0) CALL toc('Write matrix')
  
  
  CALL mpi_finalize(ierr)

END PROGRAM covmatmul
