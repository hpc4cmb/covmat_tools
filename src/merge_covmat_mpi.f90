! This program is used to combine detector covariance matrices
! into a single matrix (as inverses they are additive).
!
! February 18th 2008 Reijo Keskitalo
! 
! Revisions :
! July 25th 2008 -- can now merge matrices with different nstokes parameters
! May 11th 2010 -- simplified calling structure, can also reorder the matrix
! April 20 2011 -- Added --inverted flag to combine inverted matrices

PROGRAM merge_covmat_mpi

  USE healpix_types
  USE extension, ONLY   : nArguments, getArgument
  USE covmat_util, ONLY : tic, toc
  USE scalapack_tools

  IMPLICIT NONE

  INTEGER, PARAMETER :: nmax=100
  CHARACTER(len=filenamelen) :: argument, file_in(NMAX)
  CHARACTER(len=filenamelen) :: file_out
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: buffer_in(:,:), buffer_out(:,:)
  LOGICAL :: there, ok, inverted=.FALSE.
  INTEGER(I4B) :: npix, nstokes(100), nside, ncovmat, nstokesmax
  INTEGER(I4B) :: iargument, icovmat, nside_matrix, nstokes_matrix
  INTEGER(i4b) :: ierr, id, ntasks, ipixel, istokes, istokes2
  INTEGER(i4b) :: npix_proc, npix_proc_local, ncol_proc_local
  INTEGER(i4b) :: status(MPI_STATUS_SIZE)
  INTEGER(i4b) :: handle_in(NMAX), fileinfo_in(NMAX)
  INTEGER(i4b) :: handle_out, filemode, fileinfo_out
  INTEGER(MPI_OFFSET_KIND) :: fileoffset

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (id == 0) THEN
     ok = .TRUE.
     IF (nArguments() < 3)  THEN
        WRITE (*,'(/,a,/)') 'Usage: merge_covmat_mpi <covmat1> ' // &
             '[<covmat2> [<covmat3> ...]] [-o <covmat_out>]'
        ok = .FALSE.
     ELSE
        WRITE(*,'(/,a,i0,a,/)') &
             ' merge_covmat_mpi started with ',ntasks,' tasks.'

        file_in = ' '
        nstokes = 0
        ncovmat = 0
        file_out = 'merged_covmat.dat'
        iArgument = 1
        DO
           CALL getargument(iArgument, argument)
           IF (INDEX(argument,'-o') /= 0) THEN
              iArgument = iArgument + 1
              CALL getargument(iArgument, argument)
              file_out = TRIM(ADJUSTL(argument))
           ELSE IF (INDEX(argument,'-inv') /= 0) THEN
              inverted = .true.
              write (*,*) ' Assuming matrices are INVERTED'
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
              WRITE(*,'(a,i0,a)') &
                   'covmat', nCovmat, '  == ' // TRIM(file_in(ncovmat))
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

  CALL mpi_bcast(inverted, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(file_in, filenamelen*NMAX, MPI_CHARACTER, 0, &
       mpi_comm_world, ierr)
  CALL mpi_bcast(file_out, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, &
       ierr)
  CALL mpi_bcast(ncovmat, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)

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
  ALLOCATE(buffer_in(npix*nstokesmax, nstokesmax), &
       buffer_out(npix*nstokesmax, nstokesmax), &
       stat=ierr)
  IF (ierr /= 0) STOP 'no room for buffers'
        
  if (npix_proc_local > 0) then
     DO ipixel = npix_proc*id, npix_proc*id + npix_proc_local - 1
        buffer_out = 0.0_dp          
        ! Read and merge
        DO icovmat = 1, ncovmat
           fileoffset = int(ipixel,i8b) * npix * nstokes(icovmat)**2
           CALL mpi_file_read_at(handle_in(icovmat), fileoffset, buffer_in, &
                npix*nstokes(icovmat)**2, MPI_REAL8, status, ierr)
           
           IF (nstokes(icovmat) /= nstokesmax) THEN
              ! replace sentinel value -1 by zero for merging
              if (buffer_in(ipixel+1,1) < 0) buffer_in(ipixel+1,1) = 0.0_dp

              ! merge temperature-only matrix
              buffer_out(1:npix,1) = buffer_out(1:npix,1) &
                   + buffer_in(1:npix,1)
           ELSE
              ! replace unobserved sentinel value (-1) by zero for merging
              if (buffer_in(ipixel+1,1) < 0) then
                 do istokes = 0, nstokesmax - 1
                    buffer_in(ipixel+1+istokes*npix,istokes+1) = 0.0_dp
                 end do
              end if

              ! add to outbuffer
              DO istokes = 1, nstokesmax
                 DO istokes2 = 1, nstokesmax
                    buffer_out(1+(istokes2-1)*npix:istokes2*npix, &
                         istokes) = &
                         buffer_out(1+(istokes2-1)*npix:istokes2*npix, &
                         istokes) + &
                         buffer_in(istokes2:npix*nstokesmax:nstokesmax, &
                         istokes)
                 END DO
              END DO
           END IF
        END DO

        ! inverted covariances are not additive. Assume user wants
        ! the covariance of an average map
        if (inverted) buffer_out = buffer_out / ncovmat**2

        ! replace the missing pixels again with -1
        if (buffer_out(ipixel+1,1) == 0.0_dp) then
           do istokes = 0, nstokesmax - 1
              buffer_out(ipixel+1+istokes*npix,istokes+1) = -1.0_dp
           end do
        end if

        ! Write
        DO istokes = 0, nstokesmax-1
           fileoffset = int(ipixel + istokes*npix, i8b)*npix*nstokesmax
           CALL mpi_file_write_at(handle_out, fileoffset, &
                buffer_out(:,istokes+1), npix*nstokesmax, &
                MPI_REAL8, status, ierr)
        END DO
     END DO
  end if
  
  CALL mpi_barrier(mpi_comm_world, ierr)
  do icovmat = 1,ncovmat 
     CALL mpi_file_close(handle_in(icovmat), ierr)
  end do
  CALL mpi_file_close(handle_out, ierr)
  
  if (id == 0) call toc('merge matrix')
  
  CALL mpi_finalize(ierr)

END PROGRAM merge_covmat_mpi
