! Tool to apply a linear transformation (smoothing and/or downgrading)
! to each eigenvector separately


PROGRAM downgrade_evecs_mpi

  USE pix_tools
  USE healpix_types
  USE extension

  USE covmat_util, ONLY : tic, toc

  USE downgrade
  USE scalapack_tools

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: file_in, file_out, file_nobs, file_beam
  INTEGER(i4b) :: nside_in, nside_out, nstokes, midnside
  INTEGER(i4b) :: id, ntasks, ierr
  REAL(dp) :: FWHM, rcondlim
  INTEGER(i8b) :: npix_in, npix_out
  CHARACTER(len=filenamelen) :: argument
  LOGICAL :: use_cosine_window=.FALSE., noiseweight=.FALSE., polsmooth=.TRUE.
  LOGICAL :: smooth=.FALSE., use_beam=.false.

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)
  IF (id == 0) WRITE (*,'(/,a,i0,a,/)') &
       ' downgrade_evecs_mpi started with ',ntasks,' tasks.'
  
  IF (id == 0) CALL tic
  CALL parse_arguments
  IF (id == 0) CALL toc('Parse arguments')

  CALL get_matrix_size_dp(file_in, nside_in, nstokes)  

  IF (id == 0) WRITE (*, '(a,i8)') ' nside_in  == ', nside_in
  npix_in = nside2npix(nside_in)
  IF (id == 0) WRITE (*, '(a,i8)') '  npix_in  == ', npix_in

  IF (nside_out > 0) THEN
     npix_out = nside2npix(nside_out)
  ELSE
     nside_out = nside_in
     npix_out = npix_in
  END IF
  IF (id == 0) WRITE (*, '(a,i8)') ' nside_out == ', nside_out
  IF (id == 0) WRITE (*, '(a,i8)') '  npix_out == ', npix_out

  IF (id == 0) WRITE (*, '(a,i8)') '   nstokes == ', nstokes

  IF (id == 0) CALL tic
  CALL processEvecs()
  CALL mpi_barrier(mpi_comm_world, ierr)
  IF (id == 0) CALL toc('Process eigenvectors')

  CALL mpi_finalize(ierr)


CONTAINS


  SUBROUTINE processEvecs()

    ! Read, process and write the eigenvectors one by one
    ! mpi version reads and writes all data at once to get optimal
    ! bandwidth, may require a lot of memory

    REAL(dp), POINTER :: map_in(:, :), map_out(:, :), map_mid(:, :), &
         map_unsmoothed(:, :)
    INTEGER(i4b) :: ierr
    INTEGER(i8b) :: icol
    INTEGER(i4b) :: status(MPI_STATUS_SIZE), filemode, fileinfo_in, handle_in
    INTEGER(i4b) :: fileinfo_out, handle_out
    INTEGER(i8b) :: ncol_proc, ncol_proc_local, npix_mid
    INTEGER(MPI_OFFSET_KIND) :: fileoffset

   ! calculate the data distribution
    ncol_proc = npix_in*nstokes / ntasks
    IF (ncol_proc*ntasks < npix_in*nstokes) ncol_proc = ncol_proc + 1
    IF (id == 0) WRITE (*,'(" Processing ",i0," columns per task")') ncol_proc

    ! This task may get less pixels
    ncol_proc_local = ncol_proc
    IF (ncol_proc*(id+1) > npix_in*nstokes) THEN
       ncol_proc_local =  npix_in*nstokes - ncol_proc*id
       IF (ncol_proc_local < 1) THEN
          ncol_proc = 0
          ncol_proc_local = 0
       END IF
       WRITE (*,'(i4,": special ncol_proc_local == ",i0)') id, ncol_proc_local
    END IF

    ! allocate space

    if ( smooth ) then
       if (midnside > 0) then
          if (midnside > nside_in) stop 'midnside must be less than input nside'
          if (midnside < nside_out) stop 'midnside must be greater than or equal the output nside'
       else
          IF (nside_out < nside_in ) THEN
             midnside = nside_out * 2 
          ELSE
             midnside = nside_out
          END IF
       end if
    else
       if (midnside > 0) print *,'Warning, midnside ignored without smoothing'
       midnside = nside_out
    end if

    npix_mid = nside2npix(midnside)

    IF (ncol_proc_local > 0) THEN
       ALLOCATE(map_in(0:npix_in-1, nstokes), map_out(0:npix_out-1, nstokes), &
            map_mid(0:npix_mid-1, nstokes), map_unsmoothed(0:npix_out-1,nstokes), &
            stat=ierr)
       IF (ierr /= 0) STOP 'no room for evecs'
    ELSE
       map_in => null()
       map_mid => null()
       map_out => null()
    END IF

    ! open the mpi-io file view INPUT
    CALL mpi_info_create(fileinfo_in, ierr)
    filemode = MPI_MODE_RDONLY
    CALL mpi_file_open(mpi_comm_world, TRIM(file_in), filemode, &
         fileinfo_in, handle_in, ierr)
    IF (ierr /= 0) STOP 'Failed to open file_in'

    fileoffset = 0 ! ncol_proc * npix_in * nstokes * id
    CALL mpi_file_set_view(handle_in, fileoffset, MPI_REAL8, &
         MPI_REAL8, 'native', fileinfo_in, ierr)
    IF (ierr /= 0) STOP 'Failed to open view to file_in'

    ! open the mpi-io file view OUTPUT
    CALL mpi_info_create(fileinfo_out, ierr)
    filemode = IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE)
    CALL mpi_file_open(mpi_comm_world, TRIM(file_out), filemode, &
         fileinfo_out, handle_out, ierr)
    IF (ierr /= 0) STOP 'Failed to open file_out'

    fileoffset = 0
    CALL mpi_file_set_view(handle_out, fileoffset, MPI_REAL8, &
         MPI_REAL8, 'native', fileinfo_out, ierr)
    IF (ierr /= 0) STOP 'Failed to open view to file_out'

    ! process the eigenvectors

    DO icol = 1, ncol_proc_local
       IF (id == 0 .AND. MODULO(icol, ncol_proc/10) == 1) &
            WRITE (*,'(i4,a)') icol*100/(ncol_proc), '% done'

       ! read in the eigenvector
       fileoffset = (ncol_proc * id + icol - 1)* npix_in * nstokes
       CALL mpi_file_read_at(handle_in, fileoffset, map_in, npix_in*nstokes, &
            MPI_REAL8, status, ierr)

       if (nside_in /= midnside) then
          IF (noiseweight) then
             CALL noise_weight_downgrade(file_nobs, map_in, map_mid, nside_in, &
                  midnside, nstokes, rcondlim)
          else
             CALL average_downgrade(map_in, map_mid, nside_in, midnside, &
                  nstokes)
          end IF
       else
          map_mid(:, :) = map_in(:, :)
       end if
 
       if (smooth) then
          if (midnside /= nside_out) then
             call average_downgrade(map_mid, map_unsmoothed, midnside, &
                  nside_out, nstokes)
          else
             map_unsmoothed(:, :) = map_mid(:, :)
          end if

          IF (use_beam) THEN
             CALL windowfile_downgrade(file_beam, map_mid, map_out, midnside, &
                  nside_out, nstokes, .false., .true., .true., .false.)
          ELSE IF (use_cosine_window) THEN
             CALL apodized_downgrade(map_mid, map_out, midnside, nside_out, &
                  nstokes, .false., .true., .true., .false.)
          ELSE
             CALL smooth_downgrade(fwhm, map_mid, map_out, midnside, nside_out, &
                  nstokes, .false., .true., .true., .false.)
          END IF

          if (nstokes == 3 .and. .not. polsmooth) then
             map_out(:, 2:3) = map_unsmoothed(:, 2:3)
          end if

          where (map_unsmoothed == 0) map_out = 0
       else
          map_out(:, :) = map_mid(:, :)
       end if

       fileoffset = (ncol_proc * id + icol - 1) * npix_out * nstokes
       CALL mpi_file_write_at(handle_out, fileoffset, map_out, &
            npix_out*nstokes, MPI_REAL8, status, ierr)
       IF (ierr /= 0) STOP 'Failed to write into file_out'
    END DO

    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_file_close(handle_in, ierr)

    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_file_close(handle_out, ierr)

    if ( associated(map_in) ) deallocate(map_in, map_mid, map_out)

    IF (id == 0) WRITE (*,*) 'Done!'

  END SUBROUTINE processEvecs



  SUBROUTINE parse_arguments()

    LOGICAL :: there, ok
    integer :: iargument

    IF (id == 0) THEN
       ok = .FALSE.
       DO
          IF (nArguments() < 2) THEN
             WRITE (*,'(/,a,/)') &
                  '  Usage: downgrade_evecs_mpi <evec-file> <nside_out> ' // &
                  '[-nobs <Nobs_matrices>] [-fwhm <FWHM>] ' // &
                  '[-beam <beamfile>] [-o <out>] ' // &
                  '[-nopolsmooth] [-rcond <rcond_lim>] [-midnside <midnside>]'
             exit
          END IF

          CALL getArgument(1, file_in)
          INQUIRE(file=TRIM(file_in), exist=there)
          IF (.NOT. there) THEN
             WRITE (*,*) 'ERROR: eigenvector file does not exist'
             EXIT
          END IF

          CALL getArgument(2, argument)
          READ (argument, *, iostat=ierr) nside_out
          if (ierr /= 0) then
             print *,'Failed to parse nside_out from ' // trim(argument)
             stop
          end if
          WRITE (*,*) ' Nside out == ', nside_out

          noiseWeight = .FALSE.
          smooth = .FALSE.
          polsmooth = .TRUE.
          midnside = -1
          use_cosine_window = .FALSE.
          use_beam = .false.
          file_out = '!out.fits'
          FWHM = 0
          rcondlim = 0
          
          iargument = 3
          DO
             IF (iargument > narguments()) EXIT
             
             CALL getArgument(iargument, argument)
             IF (INDEX(argument, '-nobs') /= 0) THEN
                iArgument = iArgument + 1           
                CALL getArgument(iargument, argument)
                noiseWeight = .TRUE.
                file_nobs = TRIM(argument)
                INQUIRE(file=TRIM(file_nobs), exist=there)
                IF (.NOT. there) STOP 'file_nobs does not exist'
                WRITE (*,'(a)') ' Doing noise weighted averaging according ' &
                     // 'to ' // TRIM(file_nobs)
             ELSE IF (INDEX(argument, '-beam') /= 0) THEN
                if (smooth) stop 'Cannot apply more then one smoothing window'
                smooth = .TRUE.
                use_beam = .TRUE.
                iArgument = iArgument + 1           
                CALL getArgument(iargument, argument)
                file_beam = TRIM(argument)
                INQUIRE(file=TRIM(file_beam), exist=there)
                IF (.NOT. there) STOP 'beamfile does not exist'
                WRITE (*,'(a)') ' Doing harmonic smoothing according to ' &
                     // TRIM(file_beam)
             ELSE IF (INDEX(argument, '-fwhm') /= 0) THEN
                iArgument = iArgument + 1           
                CALL getArgument(iargument, argument)
                READ (argument,*,iostat=ierr) FWHM
                IF (ierr /= 0) STOP 'Unable to parse FWHM'
                IF (FWHM > 1E-6) THEN
                   if (smooth) stop 'Cannot apply more then one smoothing window'
                   smooth = .TRUE.
                   WRITE (*,'(a,g15.5,a)') ' Doing harmonic smoothing with ' &
                        // 'FWHM == ', FWHM, ''' beam'
                ELSE IF (FWHM < -1E-6) THEN
                   if (smooth) stop 'Cannot apply more then one smoothing window'
                   smooth          = .TRUE.
                   use_cosine_window = .TRUE.
                   WRITE (*,'(a,g15.5,a)') &
                        ' Doing harmonic smoothing with cosine window'
                ELSE
                   WRITE (*,'(a)') ' FWHM=0 :  naive pixel averaging'
                END IF
             ELSE IF (INDEX(argument, '-midnside') /= 0) THEN
                iArgument = iArgument + 1           
                CALL getArgument(iargument, argument)
                READ (argument,*,iostat=ierr) midnside
                IF (ierr /= 0) STOP 'Unable to parse midnside'
                WRITE (*,'(a,i8)') ' Using intermediate nside = ', midnside
             ELSE IF (INDEX(argument, '-rcond') /= 0) THEN
                iArgument = iArgument + 1           
                CALL getArgument(iargument, argument)
                READ (argument,*,iostat=ierr) rcondlim
                IF (ierr /= 0) STOP 'Unable to parse rcondlim'
                WRITE (*,'(a,es15.5)') ' Discarding pixels with rcond less ' &
                     // 'than ', rcondlim
             ELSE IF (INDEX(argument, '-o') /= 0) THEN
                iArgument = iArgument + 1           
                CALL getArgument(iargument, argument)
                file_out = TRIM(argument)
             ELSE IF (INDEX(argument, '-nopolsmooth') /= 0) THEN
                polsmooth = .FALSE.
                write (*,'(a)') ' ONLY smoothing the temperature component.'
             ELSE
                WRITE (*,*) 'Unrecognized argument: ' // TRIM(argument)
                exit
             END IF

             iargument = iargument + 1
          END DO

!          if (noiseweight .and. (nside_in == nside_out)) then
!             print *,'Cannot noise-weight without change in resolution'
!             exit
!          end if

          WRITE (*,'(a)') &
               ' Saving downgraded eigenvectors into '//TRIM(file_out)
          ok = .TRUE.
          EXIT
       end DO

    END IF

    CALL mpi_bcast(ok, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    IF (.NOT. ok) THEN
       CALL mpi_finalize(ierr)
       STOP
    END IF
    
    CALL mpi_bcast(file_in, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, ierr)         
    !CALL mpi_bcast(nside_in, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
    !npix_in = nside2npix(nside_in)
    CALL mpi_bcast(midnside, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(nside_out, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
    npix_out = nside2npix(nside_out)
    !CALL mpi_bcast(nstokes, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(noiseweight, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(smooth, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(use_beam, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(polsmooth, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(use_cosine_window, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(file_nobs, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, ierr)         
    CALL mpi_bcast(file_beam, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, ierr)         
    CALL mpi_bcast(FWHM, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(rcondlim, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
    CALL mpi_bcast(file_out, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, ierr)         
    
  END SUBROUTINE parse_arguments


END PROGRAM downgrade_evecs_mpi
