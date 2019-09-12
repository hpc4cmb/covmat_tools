! Add covariance matrices
!
! May 6 2013 Reijo Keskitalo
! 
! Revisions :

PROGRAM add_covmat_mpi

  USE healpix_types
  USE extension, ONLY   : nArguments, getArgument
  USE covmat_util, ONLY : tic, toc
  USE scalapack_tools

  IMPLICIT NONE

  INTEGER, PARAMETER :: nmax=100
  CHARACTER(len=filenamelen) :: argument, file_in(NMAX)
  CHARACTER(len=filenamelen) :: file_out
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), ALLOCATABLE :: buffer_in(:), buffer_out(:)
  LOGICAL :: there, ok
  INTEGER(I4B) :: npix, nstokes(100), nside, ncovmat, nstokesmax
  INTEGER(I4B) :: iargument, icovmat, nside_matrix, nstokes_matrix
  INTEGER(I4b) :: ierr, id, ntasks, ipixel, istokes, idiag
  INTEGER(I4b) :: npix_proc, npix_proc_local, ncol_proc_local
  INTEGER(I4b) :: status(MPI_STATUS_SIZE)
  INTEGER(I4b) :: handle_in(NMAX), fileinfo_in(NMAX)
  INTEGER(I4b) :: handle_out, filemode, fileinfo_out
  REAL(dp) :: factors(NMAX)
  INTEGER(MPI_OFFSET_KIND) :: fileoffset

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (id == 0) THEN
     ok = .TRUE.
     IF (nArguments() < 3)  THEN
        WRITE (*,'(/,a,/)') 'Usage: add_covmat_mpi <covmat1> <factor1> ' // &
             '[<covmat2> <factor2> [<covmat3> <factor3>...]] [-o <covmat_out>]'
        ok = .FALSE.
     ELSE
        WRITE(*,'(/,a,i0,a,/)') &
             ' add_covmat_mpi started with ',ntasks,' tasks.'

        file_in = ' '
        nstokes = 0
        ncovmat = 0
        file_out = 'summed_covmat.dat'
        iArgument = 1
        DO
           CALL getargument(iArgument, argument)
           IF (INDEX(argument,'-o') /= 0) THEN
              iArgument = iArgument + 1
              CALL getargument(iArgument, argument)
              file_out = TRIM(ADJUSTL(argument))
           ELSE
              CALL getargument(iArgument, argument)
              nCovmat = nCovmat + 1
              file_in(nCovmat) = TRIM(argument)
              INQUIRE(file=file_in(nCovmat), exist=there)
              IF (.NOT. there) THEN
                 WRITE (*,'("ERROR: ",a," does not exist!")') &
                      TRIM(file_in(ncovmat))
                 ok = .FALSE.
                 EXIT
              END IF

              iargument = iargument + 1
              call getargument(iargument, argument)
              read(argument, *, iostat=ierr) factors(ncovmat)
              if (ierr /= 0) then
                 write (*,*) 'ERROR: Failed to parse ', trim(argument), &
                      ' for scaling.'
                 ok = .FALSE.
                 EXIT
              end if

              WRITE(*,'(a,i0,a,g20.10)') &
                   'covmat', nCovmat, '  = ' // TRIM(file_in(ncovmat)) &
                   // ', factor = ',factors(ncovmat)
           END IF

           iargument = iargument + 1
           IF (iargument > narguments()) EXIT
        END DO

       IF (ok) WRITE (*,*) ' Writing to ' // TRIM(file_out)
     END IF
  END IF

  CALL mpi_bcast(ok, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
  IF (.NOT. ok) THEN
     CALL mpi_finalize(ierr)
     STOP
  END IF

  CALL mpi_bcast(ncovmat, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(file_in, filenamelen*ncovmat, MPI_CHARACTER, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(factors, ncovmat, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(file_out, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, ierr)

  do icovmat = 1, ncovmat
     CALL get_matrix_size_dp(file_in(icovmat), &
          nside_matrix, nstokes_matrix)
     IF (nside_matrix <= 0) THEN
        print *,' ERROR: Unable to determine resolution of ' &
             // trim(file_in(icovmat))
        call mpi_abort(mpi_comm_world, -1, ierr)
     END IF
     IF (icovmat == 1) THEN
        nside = nside_matrix
     ELSE IF (nside /= nside_matrix) THEN
        if (id == 0) WRITE (*,*) 'ERROR: nside of the matrices does not match'
        call mpi_abort(mpi_comm_world, -1, ierr)
     END IF

     nstokes(icovmat) = nstokes_matrix
     IF (id == 0) THEN
        WRITE (*,'(a,i0,a,i4)') 'nstokes', iCovmat, '  == ', nstokes_matrix
     ENDIF
  end do

  npix = 12*nside**2
  nstokesmax = MAXVAL(nstokes(1:ncovmat))

  call tic()

  ! determine the data distribution
  npix_proc = CEILING(DBLE(npix)/ntasks)
  IF (id == 0) &
       WRITE (*,'("Distributing matrix ",i0," pixels per task")') npix_proc
  npix_proc_local = npix_proc
  IF ((id+1)*npix_proc > npix) THEN
     npix_proc_local = npix - id*npix_proc
     IF (npix_proc_local < 0) npix_proc_local = 0
  END IF
  ncol_proc_local = npix_proc_local * nstokesmax

  ! open the mpi-io file view INPUT
  DO icovmat = 1, ncovmat 
     CALL mpi_info_create(fileinfo_in(icovmat), ierr)
     filemode = MPI_MODE_RDONLY
     CALL mpi_file_open(mpi_comm_world, file_in(icovmat), filemode, &
          fileinfo_in(icovmat), handle_in(icovmat), ierr)
     IF (ierr /= 0) STOP 'Failed to open file_in'
     
     !ncol_proc = npix_proc * nstokes(icovmat)
     !fileoffset = ncol_proc * npix * nstokes(icovmat) * id
     !if (fileoffset > (npix*nstokes(icovmat))**2-1) fileoffset = 0
     fileoffset = 0
     CALL mpi_file_set_view(handle_in(icovmat), fileoffset, MPI_REAL8, &
          MPI_REAL8, 'native', fileinfo_in(icovmat), ierr)
     IF (ierr /= 0) STOP 'Failed to open view to file_in'
  END DO

  ! open the mpi-io file view OUTPUT
  CALL mpi_info_create(fileinfo_out, ierr)
  filemode = IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE)
  CALL mpi_file_open(mpi_comm_world, file_out, filemode, &
       fileinfo_out, handle_out, ierr)
  IF (ierr /= 0) STOP 'Failed to open file_out'
        
  fileoffset = 0
  CALL mpi_file_set_view(handle_out, fileoffset, MPI_REAL8, &
       MPI_REAL8, 'native', fileinfo_out, ierr)
  IF (ierr /= 0) STOP 'Failed to open view to file_out'
        
  ! loop over local portions in all the matrices, merge and write out
  ALLOCATE(buffer_in(npix*nstokesmax), buffer_out(npix*nstokesmax), stat=ierr)
  IF (ierr /= 0) STOP 'no room for buffers'
        
  if (npix_proc_local > 0) then
     loop_pixel : DO ipixel = npix_proc*id, npix_proc*id + npix_proc_local - 1
        ! Read and merge, assume that the matrix is in block format
        loop_stokes : DO istokes = 1, nstokesmax
           idiag = ipixel + 1 + (istokes-1)*npix ! diagonal index
           buffer_out = 0.0_dp
           loop_covmat : DO icovmat = 1, ncovmat
              
              if (istokes > nstokes(icovmat)) cycle loop_covmat

              fileoffset = int(ipixel + (istokes-1)*npix, i8b)*npix*nstokes(icovmat)
              CALL mpi_file_read_at(handle_in(icovmat), fileoffset, buffer_in, &
                   npix*nstokes(icovmat), MPI_REAL8, status, ierr)
           
              ! replace unobserved sentinel value (-1) by zero for merging
              if (buffer_in(idiag) < 0) buffer_in(idiag) = 0.0_dp

              ! scale
              if (factors(icovmat) /= 1) buffer_in = buffer_in * factors(icovmat)

              IF (nstokes(icovmat) /= nstokesmax) THEN
                 ! Read nstokes==1, write nstokes==3
                 buffer_out(1:npix) = buffer_out(1:npix) + buffer_in(1:npix)
              ELSE
                 buffer_out = buffer_out + buffer_in
              END IF
           END DO loop_covmat
        
           ! replace the missing pixels again with -1
           if (buffer_out(idiag) == 0.0_dp) buffer_out(idiag) = -1.0_dp

           ! Write
           fileoffset = int(ipixel + (istokes-1)*npix, i8b)*npix*nstokesmax
           CALL mpi_file_write_at(handle_out, fileoffset, buffer_out, &
                npix*nstokesmax, MPI_REAL8, status, ierr)

        END DO loop_stokes
     END DO loop_pixel
  end if
  
  CALL mpi_barrier(mpi_comm_world, ierr)
  do icovmat = 1,ncovmat 
     CALL mpi_file_close(handle_in(icovmat), ierr)
  end do
  CALL mpi_file_close(handle_out, ierr)
  
  if (id == 0) call toc('add matrix')
  
  CALL mpi_finalize(ierr)

END PROGRAM add_covmat_mpi
