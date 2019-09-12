! This simple code adds an offset to the covmat
! diagonal.
 
PROGRAM add_to_diagonal

  USE healpix_types
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file, outfilename
  INTEGER(i4b)               :: npix, nstokes, nside, rec_len
  INTEGER(i8b)               :: icol
  INTEGER, PARAMETER         :: covmat_unit=55
  REAL(dp), POINTER          :: covmat(:)
  REAL(dp)                   :: offset, pol_offset
  LOGICAL                    :: there
  INTEGER(i4b)               :: ierr

  IF (nArguments() < 5)  THEN
     WRITE (*,*) 'Usage: add_to_diagonal <covmat> <nside> <nstokes> ' // &
          '<TT offset> <Pol offset> [<outfile>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *, iostat=ierr) nside
     if (ierr /= 0) then
        write (*,'(a)') 'Failed to parse ' // trim(argument) // ' for nside'
        stop
     end if
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(3, argument)
     READ (argument, *, iostat=ierr) nstokes
     if (ierr /= 0) then
        write (*,'(a)') 'Failed to parse ' // trim(argument) // ' for nstokes'
        stop
     end if
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes

     CALL getArgument(4, argument)
     READ (argument, *, iostat=ierr) offset
     if (ierr /= 0) then
        write (*,'(a)') 'Failed to parse ' // trim(argument) // ' for II offset'
        stop
     end if
     WRITE (*,*) ' Adding ', offset, ' to II diagonal'


     CALL getArgument(5, argument)
     READ (argument, *, iostat=ierr) pol_offset
     if (ierr /= 0) then
        write (*,'(a)') 'Failed to parse ' // trim(argument) // ' for Pol offset'
        stop
     end if
     WRITE (*,*) ' Adding ', pol_offset, ' to Pol diagonal'


     outfilename = 'offset_'//TRIM(ADJUSTL(covmat_file))
     IF (nArguments() > 5) THEN
        CALL getArgument(6, argument)
        outfilename = TRIM(ADJUSTL(argument))
     END IF
     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  npix = 12*nside**2

  ! read, add offset and write. one column at a time
  CALL tic
  ALLOCATE(covmat(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for covmat'
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)
  OPEN(unit=covmat_unit+1, file=TRIM(outfilename), status='replace', &
       form='unformatted', access='direct', recl=rec_len)

  DO icol = 0, npix - 1
     READ(covmat_unit, rec=icol+1) covmat
     covmat(icol) = covmat(icol) + offset
     WRITE(covmat_unit+1, rec=icol+1) covmat
  END DO

  DO icol = npix, npix*nstokes - 1
     READ(covmat_unit, rec=icol+1) covmat
     covmat(icol) = covmat(icol) + pol_offset
     WRITE(covmat_unit+1, rec=icol+1) covmat
  END DO

  CLOSE(covmat_unit)
  CLOSE(covmat_unit+1)
  CALL toc('add_to_diagonal')

END PROGRAM add_to_diagonal
