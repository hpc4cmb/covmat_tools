! Simple interface to PBLAS routine PDGEMM to
! multiply two block cyclic distributed double precision matrices
!
! 2008-10-26 Reijo Keskitalo


PROGRAM evecs2covmat_mpi

  USE healpix_types
  USE pix_tools, ONLY : nside2npix, npix2nside
  USE extension, ONLY : getArgument, nArguments
  USE scalapack_tools

  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  TYPE(scalapack_env) :: spenv
  TYPE(scalapack_matrix_dp) :: evecs, covmat

  CHARACTER(len=filenamelen) :: argument, file_evecs, file_evals, file_out
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: evals(:)
  LOGICAL(lgt) :: there, invert=.false., ok, sqroot=.false., notortho=.false.
  INTEGER(i4b) :: ierr, id, ntasks
  INTEGER(i4b) :: ielement, blocksize
  INTEGER(i4b) :: nstokes, nside, nside2, iargument
  REAL(dp) :: limit=1E-30
  INTEGER(i8b) :: nelem, irow, icol, npix_tot_in, npix_tot_out, npix, &
       icolumn_tot

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)
  IF (id == 0) WRITE (*,'(/,a,i0,a,/)') &
       ' evecs2covmat_mpi started with ',ntasks,' tasks.'
  
  CALL parse_arguments

  CALL get_matrix_size_dp(file_evecs, nside, nstokes, npixtot=npix, nelem=nelem)
  if (id == 0) write (*,*) 'Evec file has ',nelem,' elements'

  if (npix_tot_out > 0) then
     npix_tot_in = nelem / npix_tot_out
  else
     npix_tot_in = npix
     npix_tot_out = npix
  end if

  IF (id == 0) WRITE (*,'(" npix_tot_in  == ",i0)') npix_tot_in
  IF (id == 0) WRITE (*,'(" npix_tot_out == ",i0)') npix_tot_out

  IF (id == 0) THEN
     WRITE (*,'("Block == ",i4)') blocksize
     WRITE (*,*) ' writing to ' // TRIM(file_out)
  END IF

  CALL init_scalapack(spenv, blocksize)

  CALL mpi_barrier(mpi_comm_world, ierr)

  ! read in the matrices
  ALLOCATE(evals(npix_tot_in), stat=ierr)
  IF (ierr /= 0) STOP 'no room for covmat'

  IF (ID == 0) CALL tic
     
  IF (id == 0) THEN
     OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(file_evals)), status='old',&
          form='formatted')
     READ(covmat_unit, *, iostat=ierr) evals
     IF (ierr /= 0) STOP 'Failed to read eigenvalues from disk'
     CLOSE(covmat_unit)

     WRITE (*,'(/,a)') 'Eigenvalues :'
     WRITE (*,'(es20.10)') evals(1:10)
     WRITE (*,*) '...'
     WRITE (*,'(es20.10)') evals(npix_tot_in-10:npix_tot_in-1)
     WRITE (*,'(/)')
  END IF
  call mpi_bcast(evals, npix_tot_in, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)

  limit = MAXVAL(evals) * limit
  WHERE(evals < limit) evals = 0
  IF (invert) THEN
     IF (id == 0) WRITE (*,*) 'Inverting the eigenvalues'
     where(evals > 0) evals = 1.0_dp / evals
  END IF

  if (sqroot) evals = sqrt(evals)

  CALL mpi_barrier(mpi_comm_world, ierr)

  ! npix_tot_out and npix_tot_in were swapped on April 6 2012 to
  ! fix a bug in handling of non-square eigenvector files.
  CALL read_scalapack_matrix_dp(spenv, evecs, file_evecs, npix_tot_out, npix_tot_in)

  CALL mpi_barrier(mpi_comm_world, ierr)

  IF (ID == 0) CALL toc('Read matrices')
  IF (ID == 0) CALL tic
  
  if (notortho) then
     ! Eigenvectors are not assumed to be orthogonal, so the recomposition
     ! must be computed one vector at a time.
     call initialize_scalapack_matrix_dp( &
          spenv, covmat, int(npix_tot_out, i4b), int(npix_tot_out, i4b))
     do icol = 1,npix_tot_in
        call pdger(npix_tot_out, npix_tot_out, evals(icol), &
             evecs%data, 1, icol, evecs%desc, 1, &
             evecs%data, 1, icol, evecs%desc, 1, &
             covmat%data, 1, 1, covmat%desc)
     end do
  else
     ! apply the eigenvalues by looping over the local array
     ! and converting locations into global column numbers
     ielement = 0
     evals = SQRT(evals)
     DO icol = 0, evecs%ncol - 1
        DO irow = 0, evecs%nrow - 1
           ielement = ielement + 1        
           icolumn_tot = &
                (INT(icol/spenv%block, i8b)*spenv%npcol &
                + spenv%mycol)*spenv%block &
                + MODULO(icol, int(spenv%block, i8b))
           
           evecs%data(ielement) = evecs%data(ielement) * evals(icolumn_tot+1)
        END DO
     END DO
  
     IF (ID == 0) CALL toc('Apply eigenvalues')
     IF (ID == 0) CALL tic
     
     ! Multiply
     CALL multiply_scalapack_matrix_dp(spenv, evecs, evecs, covmat, &
          (/.FALSE., .TRUE./))
  end if

  IF (ID == 0) CALL toc('Multiplication')     
  IF (ID == 0) CALL tic

  IF (covmat%nrow_tot /= npix_tot_out .or. covmat%ncol_tot /= npix_tot_out) then
     write (*,*) 'ERROR: Multiplied matrix is ', &
          covmat%nrow_tot, ' x ', covmat%ncol_tot
  end IF
  
  CALL write_scalapack_matrix_dp(spenv, covmat, file_out)
  if (id == 0) write (*,'(" Recomposed matrix saved in ",a)') trim(file_out)

  IF (ID == 0) CALL toc('Write matrix')     

  CALL mpi_finalize(ierr)


CONTAINS


  SUBROUTINE parse_arguments

    ok = .TRUE.

    IF (narguments() < 2)  THEN
       IF (id == 0) WRITE(*,*) &
            'Usage: evecs2covmat_mpi ' // &
            '<evals> <evecs> [-npix npix_tot_out] ' // &
            '[--sqrt] [--inv] [--notortho] [-o <outfile>] ' // &
            '[-lim <limit>] [-b <blocksize>]'
       ok = .FALSE.
    ELSE IF (id == 0) THEN
       CALL getargument(1, argument)
       file_evals = TRIM(argument)
       INQUIRE(file=TRIM(file_evals), exist=there)
       IF (.NOT. there) THEN
          write (*,*) 'file_evals does not exist!'
          ok = .false.
       else
          WRITE (*,'(a)') ' Reading eigenvalues from ' // TRIM(file_evals)
       end IF

       CALL getargument(2, argument)
       file_evecs = TRIM(argument)
       INQUIRE(file=TRIM(file_evecs), exist=there)
       IF (.NOT. there) then
          write (*,*) 'file_evecs does not exist!'
          ok = .false.
       else
          WRITE (*,'(a)') ' Reading eigenvectors from ' // TRIM(file_evecs)
       end IF
          
       file_out = 'out.dat'
       blocksize = 32
       invert = .FALSE.
       npix_tot_out = -1
       iArgument   = 3
       DO
          IF (iargument > narguments()) EXIT

          CALL getargument(iargument, argument)
          IF (INDEX(argument, '-o') /= 0) THEN
             iargument = iargument + 1
             CALL getargument(iargument, argument)
             file_out = TRIM(ADJUSTL(argument))
          ELSE IF (INDEX(argument, '-lim') /= 0) THEN
             iargument = iargument + 1
             CALL getargument(iargument, argument)
             READ(argument, *, iostat=ierr) limit
             IF (ierr /= 0) THEN
                write (*,*) 'Failed to parse limit from "' &
                     // trim(argument) // '"'
                ok = .FALSE.
                EXIT
             END IF
             WRITE (*,*) 'Using threshold eval > evalmax * ',limit
          ELSE IF (INDEX(argument, '-npix') /= 0) THEN
             iargument = iargument + 1
             CALL getargument(iargument, argument)
             READ(argument, *, iostat=ierr) npix_tot_out
             IF (ierr /= 0) THEN
                write (*,*) 'Failed to parse npix from "' &
                     // trim(argument) // '"'
                ok = .FALSE.
                EXIT
             END IF
             WRITE (*,*) 'npix_tot_out == ', npix_tot_out
          ELSE IF (INDEX(argument, '-b') /= 0) THEN
             iargument = iargument + 1
             CALL getargument(iargument, argument)
             READ(argument, *, iostat=ierr) blocksize
             IF (ierr /= 0) THEN
                write (*,*) 'Failed to parse blocksize from "' &
                     // trim(argument) // '"'
                ok = .FALSE.
                EXIT
             END IF
          ELSE IF (INDEX(argument, '-inv') /= 0) THEN
             invert = .TRUE.
             WRITE (*,*) 'Inverting while composing'
          ELSE IF (INDEX(argument, '-sqrt') /= 0) THEN
             sqroot = .TRUE.
             WRITE (*,*) 'Composing the square root'
          ELSE IF (INDEX(argument, '-notortho') /= 0) THEN
             notortho = .TRUE.
             WRITE (*,*) 'Will not assume eigenvectors to be orthogonal'
          ELSE
             WRITE (*,*) 'Unrecognized command line option: ', TRIM(argument)
             ok = .false.
          END IF

          iargument = iargument + 1
       END DO
    END IF


    CALL mpi_bcast(ok, 1, mpi_logical, 0, mpi_comm_world, ierr)
    IF (.NOT. ok) THEN
       CALL mpi_barrier(mpi_comm_world, ierr)
       CALL mpi_finalize(ierr)
       STOP
    END IF

    CALL mpi_bcast(sqroot, 1, mpi_logical, 0 ,mpi_comm_world, ierr)
    CALL mpi_bcast(notortho, 1, mpi_logical, 0 ,mpi_comm_world, ierr)
    CALL mpi_bcast(file_evals,  filenamelen,mpi_character,0,mpi_comm_world,ierr)
    CALL mpi_bcast(file_evecs,  filenamelen,mpi_character,0,mpi_comm_world,ierr)
    CALL mpi_bcast(nside2,      1,          mpi_integer,  0,mpi_comm_world,ierr)
    CALL mpi_bcast(nstokes,     1,          mpi_integer,  0,mpi_comm_world,ierr)
    CALL mpi_bcast(npix_tot_out,1,          mpi_integer,  0,mpi_comm_world,ierr)
    CALL mpi_bcast(blocksize,   1,          mpi_integer,  0,mpi_comm_world,ierr)
    CALL mpi_bcast(limit,       1,  mpi_double_precision, 0,mpi_comm_world,ierr)
    CALL mpi_bcast(invert,      1,          mpi_logical,  0,mpi_comm_world,ierr)
    CALL mpi_bcast(file_out,    filenamelen,mpi_character,0,mpi_comm_world,ierr)

  END SUBROUTINE parse_arguments


END PROGRAM evecs2covmat_mpi
