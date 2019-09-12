! This code reads in the eigenvectors and eigenvalues of a
! madam pixel-pixel covariance matrix.
! It then outputs a simulated noise map conforming to
! the eigenvectors and eigenvalues.
!
! April 28, 2011 Reijo Keskitalo

PROGRAM generate_noisemap

  USE healpix_types
  USE fitstools, ONLY   : output_map, getsize_fits
  USE head_fits, ONLY   : add_card, write_minimal_header
  USE pix_tools, ONLY   : convert_ring2nest
  USE rngmod
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments
  USE scalapack_tools

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, evec_file, eval_file
  CHARACTER(len=filenamelen) :: mapfile, outfile
  CHARACTER(len=80)          :: header(1000)
  INTEGER(i4b)               :: npix, nstokes, nside
  INTEGER(i4b)               :: seed, ieval, ievec, ireal
  INTEGER, PARAMETER         :: covmat_unit=55
  REAL(dp), POINTER          :: evecs(:,:), evals(:), map(:,:), mapvec(:)
  REAL(dp)                   :: limit=1E20
  LOGICAL                    :: there, invert=.FALSE., ok
  INTEGER(i4b)               :: ierr, nhit, iargument
  TYPE(planck_rng)           :: rnghandle
  REAL(dp)                   :: gauss

  INTEGER(i4b) :: npix_proc, npix_proc_local, id, nreal, ntasks
  INTEGER(i4b) :: status(MPI_STATUS_SIZE)
  INTEGER(i4b) :: filehandle, fileinfo, filemode
  INTEGER(MPI_OFFSET_KIND) :: fileoffset

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (id == 0) WRITE (*,'(a,i0,a)') 'generate_noisemap_mpi started with ', &
       ntasks, ' tasks'

  IF (id == 0) THEN
     ok = .TRUE.
     IF (nArguments() < 2)  THEN
        WRITE (*,'(/,a,/)') ' Usage: generate_noisemap_mpi ' &
             // '<covmat-evecs> <covmat-evals> ' &
             // '[-o <map>] [-lim <limit>] [-inv] [-n <iterations>]'
        ok = .FALSE.
     ELSE
        CALL getArgument(1, argument)
        evec_file = TRIM(argument)
        INQUIRE(file=TRIM(evec_file), exist=there)
        IF (.NOT. there) THEN
           WRITE (*,*) TRIM(evec_file) // ' does not exist!'
           ok = .FALSE.
        ELSE
           WRITE (*,'(a)') ' evecs == ' // TRIM(evec_file)
        END IF

        CALL getArgument(2, argument)
        eval_file = TRIM(argument)
        INQUIRE(file=TRIM(eval_file), exist=there)
        IF (.NOT. there) THEN
           WRITE (*,*) TRIM(eval_file) // ' does not exist!'
           ok = .FALSE.
        ELSE
           WRITE (*,'(a)') ' evals == ' // TRIM(evec_file)
        END IF
        
        iargument = 3
        limit = 0.0
        invert = .FALSE.
        mapfile = '!noisemap.fits'
        nreal = 1
        DO
           IF (nArguments() < iargument) EXIT
           CALL getArgument(iargument, argument)
           
           IF (INDEX(argument, '-o') /= 0) THEN
              iArgument = iArgument + 1
              CALL getArgument(iArgument, argument)
              mapfile = TRIM(ADJUSTL(argument))
           ELSE IF (INDEX(argument, '-lim') /= 0) THEN
              iArgument = iArgument + 1
              IF (nArguments() < iargument) &
                   STOP 'Must supply argument with -lim'
              CALL getArgument(iArgument, argument)
              READ(argument, *, iostat=ierr) limit
              IF (ierr /= 0) STOP 'Cannot parse limit'
              WRITE (*,*) 'Using threshold eval > evalmax * ',limit
           ELSE IF (INDEX(argument, '-n') /= 0) THEN
              iArgument = iArgument + 1
              IF (nArguments() < iargument) &
                   STOP 'Must supply argument with -n'
              CALL getArgument(iArgument, argument)
              READ(argument, *, iostat=ierr) nreal
              IF (ierr /= 0) STOP 'Cannot parse number of iterations'
              WRITE (*,'(a,i0,a)') ' Generating ',nreal,' noise maps'
           ELSE IF (INDEX(argument, '-inv') /= 0) THEN
              invert = .TRUE.
              WRITE (*,*) 'Inverting eigenvalues'
           ELSE
              WRITE (*,*) 'Unrecognized command line option: ',TRIM(argument)
              STOP
           END IF

           iargument = iargument + 1
        END DO

        WRITE (*,'(a)') ' outroot == ' // TRIM(mapfile)
     END IF
  END IF

  CALL mpi_bcast(ok, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
  IF (.NOT. ok) THEN
     CALL mpi_finalize(ierr)
     STOP
  END IF

  CALL mpi_bcast(evec_file, filenamelen, MPI_CHARACTER, 0, &
       mpi_comm_world, ierr)
  CALL mpi_bcast(eval_file, filenamelen, MPI_CHARACTER, 0, &
       mpi_comm_world, ierr)
  CALL mpi_bcast(limit, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(nreal, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(invert, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(mapfile, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, ierr)


  ! determine the size of the matrix and allocate space accordingly
  CALL get_matrix_size_dp(evec_file, nside, nstokes)
  IF (nside <= 0) THEN
     IF (id == 0) PRINT *,' ERROR: Unable to determine resolution of ' &
          // TRIM(evec_file)
     CALL mpi_finalize(ierr)
     STOP
  END IF

  npix = 12 * nside**2
  ALLOCATE(map(0:npix-1,nstokes), mapvec(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for maps'


  ! initialize and read in the local share of the eigenvector matrix

  npix_proc = CEILING(DBLE(npix)/ntasks)
  IF (id == 0) &
       WRITE (*,'("Distributing matrix ",i0," pixels per task")') npix_proc
  npix_proc_local = npix_proc
  IF ((id+1)*npix_proc > npix) THEN
     npix_proc_local = npix - id*npix_proc
     IF (npix_proc_local < 0) npix_proc_local = 0
  END IF

  ! open the mpi-io file view INPUT

  CALL tic

  CALL mpi_info_create(fileinfo, ierr)
  filemode = MPI_MODE_RDONLY
  CALL mpi_file_open(mpi_comm_world, evec_file, filemode, fileinfo, &
       filehandle, ierr)          
  IF (ierr /= 0) STOP 'Failed to open evec_file'
     
  fileoffset = 0
  CALL mpi_file_set_view(filehandle, fileoffset, MPI_REAL8, MPI_REAL8, &
       'native', fileinfo, ierr)
  IF (ierr /= 0) STOP 'Failed to open view to evec_file'

  ALLOCATE(evecs(npix*nstokes, npix_proc_local), stat=ierr)
  IF (ierr /= 0) STOP 'No room for evecs'

  fileoffset = INT(npix_proc,i8b) * npix * nstokes * id
  CALL mpi_file_read_at(filehandle, fileoffset, evecs, &
       npix*nstokes*npix_proc_local, MPI_REAL8, status, ierr)
  IF (ierr /= 0) STOP 'Failed to read in eigenvectors.'

  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL mpi_file_close(filehandle, ierr)
  
  IF (id == 0) CALL toc('Read eigenvectors')

  ! get the eigenvalues

  CALL tic

  ALLOCATE(evals(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for evecs'

  OPEN(unit=covmat_unit, file=TRIM(eval_file), status='old', &
       form='formatted')
  READ(covmat_unit,*) evals
  CLOSE(covmat_unit)

  WHERE(evals < MAXVAL(evals)*limit) evals = 0.0_dp
  IF (invert) THEN
     ! assume that evals correspond to N^-1
     WHERE(evals /= 0.0_dp) evals = 1.0/evals
  END IF

  IF (id == 0) CALL toc('Read eigenvalues')

  ! ready to simulate
  
  CALL tic

  ! every task generates the same random number sequence so results
  ! do not change with the number of tasks
  IF (id == 0) CALL SYSTEM_CLOCK(count=seed)
  CALL mpi_bcast(seed, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
  CALL rand_init(rnghandle, seed)
  
  loop_realizations : DO ireal = 0,nreal-1
     IF (nreal > 1) THEN
        outfile = update_filename(mapfile, ireal)
     ELSE
        outfile = mapfile
     END IF

     ! see if the file already exists and user does not want to
     ! overwrite
     IF (.not. outfile(1:1) == '!') THEN
        INQUIRE(file=TRIM(outfile), exist=there)
        IF (there) THEN
           IF (id == 0) WRITE (*,*) TRIM(outfile) // ' exists, skipping ... '
           CYCLE loop_realizations
        END IF
     END IF

     mapvec = 0.0
     ievec = 0
     DO ieval = 0, npix*nstokes - 1
        IF (evals(npix*nstokes-1) > 0.0_dp) THEN
           gauss = rand_gauss(rnghandle)
           nhit  = nhit + 1
           IF (ieval >= npix_proc*nstokes*id &
                .AND. ievec < npix_proc_local) THEN
              ! this process has the required eigenvector stored
              ievec = ievec + 1
              mapvec  = mapvec + gauss * evals(ieval)**0.5 * evecs(:,ievec)
           END IF
        END IF
     END DO
     
     CALL mpi_reduce(mapvec, map, npix*nstokes, MPI_DOUBLE_PRECISION, &
          MPI_SUM, 0, MPI_COMM_WORLD, ierr)     
     IF (ierr /= 0) STOP 'Failed to sum the noise map'
  
     IF (id == 0) CALL toc('Simulate map')

     IF (id == 0) THEN
        IF (nhit < npix*nstokes) &
             WRITE (*,*) 'Excluded ', npix*nstokes-nhit, &
             ' eigenvectors.'

        CALL tic
        
        header = ' '
        CALL write_minimal_header(header, 'MAP', nside=nside, &
             ordering='NESTED', creator='generate_noisemap', &
             randseed=seed, polar=(nstokes>1))
        
        WRITE (*,'(a)') 'Writing to ' // TRIM(outfile)

        CALL output_map(map, header, outfile)

        CALL toc('Write map')
     END IF
  END DO loop_realizations
  
  DEALLOCATE(map, mapvec, evecs, evals)

  CALL mpi_finalize(ierr)



CONTAINS



  FUNCTION update_filename(name, id)
    !
    ! add the realization number to the end, just before suffix
    !
    CHARACTER(len=filenamelen) :: update_filename
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER, INTENT(in) :: id

    INTEGER :: i

    i = SCAN(name, '.', .TRUE.)
    IF (i < 1) i = LEN_TRIM(name) + 1

    WRITE(update_filename,'(a,"_",i4.4,a)') name(1:i-1), id, TRIM(name(i:))

  END FUNCTION update_filename



END PROGRAM generate_noisemap
