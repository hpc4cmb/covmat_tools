! This program is used to combine detector covariance matrices
! into a single matrix (as inverses they are additive).
!
! February 18th 2008 Reijo Keskitalo
! 
! Revisions :
! July 25th 2008 -- can now merge matrices with different nstokes parameters

PROGRAM merge_covmat_sparse_mpi

  USE healpix_types
  USE extension, ONLY   : nArguments, getArgument
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(len=filenamelen) :: argument, covmat_file(100), pixellistfile(100)
  CHARACTER(len=filenamelen) :: outfilename
  integer(I4B) :: npix(0:100), nstokes(100)
  INTEGER(I4B) :: nside, ncovmat, nstokesmax
  INTEGER(I4B) :: iArgument, iCovmat, nhit, ipixel, icol_in, ipix
  INTEGER(I4B) :: npixmax, ind, ind2
  INTEGER(I4B), pointer :: pixellist(:, :), compressed2pix(:)
  INTEGER(I4B), pointer :: pix2compressed(:, :)
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmatIn(:), covmatOut(:)
  LOGICAL :: there, ok
  INTEGER :: iErr, iCol, id, ntasks
  INTEGER :: status(MPI_STATUS_SIZE)
  logical, pointer :: is_hit(:)
  INTEGER :: filemode, fileinfo, infilehandle(100)
  INTEGER :: outfilehandle
  INTEGER(MPI_OFFSET_KIND)   :: fileoffset

  CALL mpi_init(ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
  CALL mpi_comm_rank(mpi_comm_world, id, ierr)

  IF (id == 0) THEN
     ok = .TRUE.
     IF (nArguments() < 5)  THEN
        WRITE (*,'(/,a,/)') &
             'Usage: merge_covmat_mpi_sparse '//&
             '<pixellist1> <nstokes1> <covmat1> '// &
             '[<pixellist2> <nstokes2> <covmat2> '//&
             '[<pixellist3> <nstokes3> <covmat3> ...]] '// &
             '[-o <covmat_out>]'
        ok = .FALSE.
     ELSE
        WRITE(*,'(/,a,i0,a,/)') &
             ' merge_covmat_mpi_sparse started with ',ntasks,' tasks.'
        IF (ntasks < 2) THEN
           WRITE(*,*) ' ERROR, need as least 2 tasks'
           ok = .FALSE.
        END IF

        covmat_file = ' '
        pixellistfile = ' '
        nStokes = 0
        nCovmat = 0
        outfilename = 'merged_covmat.dat'
        iArgument = 1
        DO
           CALL getArgument(iArgument, argument)

           IF (INDEX(argument,'-o') /= 0) THEN
              CALL getArgument(iArgument+1, argument)
              outfilename = TRIM(ADJUSTL(argument))
              EXIT
           ENDIF

           nCovmat = nCovmat + 1

           pixellistfile(nCovmat) = TRIM(argument)
           INQUIRE(file=TRIM(pixellistfile(nCovmat)), exist=there)
           IF (.NOT. there) then
              write (*,*) 'pixellistfile does not exist : '&
                   //trim(pixellistfile(ncovmat))
              ok = .false.
              exit
           end IF
           WRITE(*,'(a,i0,a)') &
                'pixellist', nCovmat, '  == ' // TRIM(pixellistfile(nCovmat))

           CALL getArgument(iArgument+1, argument)
           READ(argument, *) nStokes(nCovmat)
           WRITE(*,'(a,i0,a,i4)') 'nStokes', nCovmat, ' == ', nStokes(nCovmat)

           CALL getArgument(iArgument+2, argument)
           covmat_file(nCovmat) = TRIM(argument)
           INQUIRE(file=TRIM(covmat_file(nCovmat)), exist=there)
           IF (.NOT. there) STOP 'covmat_file does not exist!'
           WRITE(*,'(a,i0,a)') &
                'covmat', nCovmat, '  == ' // TRIM(covmat_file(nCovmat))

           iArgument = iArgument + 3
           if (iArgument > nArguments()) exit
        END DO

        if (ok) WRITE (*,*) ' Writing to ' // TRIM(outfilename)
     END IF
  END IF

  CALL mpi_bcast(ok, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr)
  IF (.NOT. ok) THEN
     CALL mpi_finalize(ierr)
     STOP
  END IF

  CALL mpi_bcast(covmat_file, filenamelen*100, MPI_CHARACTER, 0, &
       mpi_comm_world, ierr)
  CALL mpi_bcast(pixellistfile, filenamelen*100, MPI_CHARACTER, 0, &
       mpi_comm_world, ierr)
  CALL mpi_bcast(outfilename, filenamelen, MPI_CHARACTER, 0, mpi_comm_world, &
       ierr)
  CALL mpi_bcast(nside,   1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(nCovmat, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(nstokes, 100, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierr)

  !nPix = 12*nSide**2
  nStokesMax = MAXVAL(nStokes)

  ! read in the pixellists
  do icovmat = 1, ncovmat
     open(file=pixellistfile(icovmat),unit=55,status='old')
     read(55,*) nside, npix(icovmat)
     if (id==0) write (*,'(a,i0,a,i0)') ' npix',icovmat,' == ',npix(icovmat)

     if (icovmat == 1) then
        ! initialize
        allocate(pixellist(12*nside**2*nstokesmax,ncovmat), stat=ierr)
        if (ierr /= 0) stop 'No room for pixellist'
        pixellist = -1
     end if

     read (55,*) pixellist(1:npix(icovmat),icovmat)
     close(55)

     if (nstokes(icovmat) == 3) then
        ! complete the mapping to include polarized pixels
        pixellist(npix(icovmat)+1:2*npix(icovmat),icovmat) = &
             pixellist(1:npix(icovmat),icovmat) + 12*nside**2
        pixellist(2*npix(icovmat)+1:3*npix(icovmat),icovmat) = &
             pixellist(1:npix(icovmat),icovmat) + 2*12*nside**2
     end if
  end do

  ! create a pixellist for the combined matrix and establish two mappings:
  ! pix2compressed -- pixel numbers to compressed pixel numbers
  ! compressed2pix -- reverse, only needed for the combined matrix
  allocate(pix2compressed(0:nstokesmax*12*nside**2-1,0:ncovmat),stat=ierr)
  if (ierr/=0) stop 'no room for pixel mapping array'
  pix2compressed = -1
  do icovmat = 1, ncovmat
     do ipixel = 1, npix(icovmat)*nstokes(icovmat)
        pix2compressed(pixellist(ipixel,icovmat),icovmat) = ipixel
     end do
  end do

  allocate(is_hit(0:12*nside**2-1),stat=ierr)
  if (ierr/=0) stop 'no room for hittable'
  is_hit = .false.
  do icovmat = 1, ncovmat
     is_hit(pixellist(1:npix(icovmat),icovmat)) = .true.
  end do
  npix(0) = count(is_hit)

  if (id == 0) then
     open(file='pixellist.dat',status='replace',unit=55)
     write (55,*) nside,npix(0)
  end if

  allocate(compressed2pix(npix(0)*nstokesMax),stat=ierr)
  if(ierr/=0) stop 'no room for compressed2pix'
  nhit = 0
  do ipix = 0, 12*nside**2-1
     if (is_hit(ipix)) then
        if (id == 0) write(55,'(i0)') ipix
        nhit = nhit + 1
        compressed2pix(nhit) = ipix
        pix2compressed(ipix,0) = nhit

        if (nstokesMax == 3) then
           compressed2pix(nhit+npix(0))   = ipix+12*nside**2
           compressed2pix(nhit+2*npix(0)) = ipix+2*12*nside**2

           pix2compressed(ipix+12*nside**2,0)   = nhit+npix(0)
           pix2compressed(ipix+2*12*nside**2,0) = nhit+2*npix(0)
        end if
     end if
  end do

  if (id == 0) then
     close(55)
     write (*,'(a,i0)') 'Merged matrix has ',npix(0),' pixels'
  end if

  npixmax = maxval(npix(1:))


  ! open the input covmat files
  CALL mpi_info_create(fileinfo, ierr)
  fileoffset = 0
  do icovmat = 1, ncovmat
     CALL mpi_file_open(mpi_comm_world,covmat_file(icovmat),MPI_MODE_RDONLY,&
          fileinfo,infilehandle(icovmat),ierr)
     CALL mpi_file_set_view(infilehandle(icovmat),fileoffset,MPI_REAL8,&
          MPI_REAL8,'native',fileinfo,ierr)
  end do

  ! open the output file
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  CALL mpi_file_open(mpi_comm_world,outfilename,filemode,fileinfo,&
       outfilehandle,ierr)
  CALL mpi_file_set_view(outfilehandle,fileoffset,MPI_REAL8,MPI_REAL8, &
       'native',fileinfo,ierr)

  ! read the matrices and write the sum
  ALLOCATE(covmatIn(npixmax*nstokesmax),covmatOut(npix(0)*nstokesmax),&
       stat=iErr)
  IF (iErr /= 0) STOP 'no room for covmat'

  DO iCol = 1, nPix(0)*nStokesMax
     IF (id == MODULO(iCol-1, ntasks)) THEN
        covmatOut = 0.0
        DO iCovmat = 1, nCovmat
           icol_in = pix2compressed(compressed2pix(icol),icovmat)
           if (icol /= icol_in) &
                write (*,'(i0,a,i0,a,i0,a,i0)') &
                id,' : icovmat=',icovmat,',icol=',icol,',icol_in=',icol_in
           IF (iCol_in > -1) THEN
              fileoffset = (iCol_in-1)*npix(icovmat)*nstokes(icovmat)
              call mpi_file_read_at(infilehandle(icovmat), &
                   fileoffset, covmatIn, &
                   npix(icovmat)*nstokes(icovmat), MPI_REAL8, status, ierr)

              ! add the read covariance matrix
              do ind = 1, npix(icovmat)*nstokes(icovmat)
                 ipixel = pixellist(ind,icovmat)
                 ind2   = pix2compressed(ipixel,0)
                 if (ind /= ind2) &
                      write (*,'(i0,a,i0,a,i0,a,i0)') &
                      id,' : ind=',ind,',ipixel=',ipixel,',ind2=',ind2
                 if (ind2 > -1) then
                    covmatOut(ind2) = covmatOut(ind2) + covmatIn(ind)
                 end if
              end do

              !covmatOut(pix2compressed(pixellist(&
              !     1:npix(icovmat)*nstokes(icovmat),icovmat),0)) = &
              !     covmatOut(pix2compressed(pixellist(&
              !     1:npix(icovmat)*nstokes(icovmat),icovmat),0)) &
              !     + covmatIn(npix(icovmat)*nstokes(icovmat))
           END IF
        END DO
        fileoffset = (iCol-1)*npix(0)*nstokesmax
        call mpi_file_write_at(outfilehandle, fileoffset, &
             covmatOut, npix(0)*nstokesmax, MPI_REAL8, status, ierr)
     END IF
  END DO

  do icovmat = 1, ncovmat
     call mpi_file_close(infilehandle(icovmat),ierr)
  end do
     call mpi_file_close(outfilehandle,ierr)

  IF (id == 0) CALL toc('merge covmat')

  DEALLOCATE(covmatIn, covmatOut)

  CALL mpi_finalize(ierr)

END PROGRAM merge_covmat_sparse_mpi
