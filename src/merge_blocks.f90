! Command line utility to merge the blocks of the covariance matrix
! after the massive parallel run.
!
! input  : fortran binary blocks from madam_NCVM run
! output : a single monolithic inverse covariance matrix file, C binary


PROGRAM merge_blocks

  USE healpix_types
  USE extension, ONLY   : getArgument, nArguments
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=80)      :: argument, covmat_file, outfilename
  INTEGER(dp)            :: npix, nstokes, nside
  INTEGER, PARAMETER     :: covmat_unit=55
  REAL(dp), POINTER      :: covmat(:)
  LOGICAL                :: there
  INTEGER(i4b)           :: ierr, icol, rec_len

  IF (nArguments() /= 3)  THEN
     STOP 'Usage: merge_blocks <first_block> <nside> <nstokes>'
  ELSE
     CALL getArgument(1, argument)
     covmat_file = TRIM(argument)
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(covmat_file)

     CALL getArgument(2, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes
  END IF

  CALL get_outfilename(covmat_file, outfilename)
  WRITE (*,*) 'Writing to ' // TRIM(outfilename)

  npix = 12*nside**2
  ALLOCATE(covmat(0:npix*nstokes-1), stat=ierr) ! need only one column
  IF (ierr /= 0) STOP 'no room for covmat'

  CALL tic
  inquire(iolength=rec_len) covmat
  OPEN(unit=covmat_unit+1, file=TRIM(ADJUSTL(outfilename)), status='replace',&
       form='unformatted', access='direct', recl=rec_len)
  icol = 0
  DO
     ! Keep reading blocks column by column
     ! and writing them to outfile
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (.NOT. there) EXIT
     OPEN(unit=covmat_unit,file=TRIM(covmat_file),status='old',&
          form='unformatted')
     IF (MODULO(icol, 100) == 0) THEN
        WRITE (*,'(/,a)') ' Reading ' // TRIM(covmat_file)
     END IF
     DO
        write (*,*) icol
        READ (covmat_unit, END=999) covmat
        icol = icol+1
        WRITE (covmat_unit+1, rec=icol) covmat
     END DO
999  CLOSE(unit=covmat_unit)
     CALL update_covmat_file(covmat_file)
     IF (MODULO(icol, 100) == 0) CALL toc('This much')
  END DO
  CLOSE(unit=covmat_unit+1)
  CALL toc('merge blocks')

  IF (icol /= npix*nstokes) write (*,'(a,i8,a,i8)') &
       ' ERROR: found ', icol, ' columns, expected ', npix*nstokes


CONTAINS


  SUBROUTINE get_outfilename(covmat_file, outfilename)
    CHARACTER(len=80) :: covmat_file, outfilename
    INTEGER           :: istart, iend

    istart = INDEX(covmat_file, '_block_')
    IF (istart < 0) STOP 'bad block number'
    iend = istart + 7 + INDEX(covmat_file(istart+7:80), '_')

    outfilename = covmat_file(1:istart) // covmat_file(iend:80)
        
  END SUBROUTINE get_outfilename


  SUBROUTINE update_covmat_file(covmat_file)
    CHARACTER(len=80) :: covmat_file, stemp
    INTEGER           :: istart, iend, blok, ilength

    istart = INDEX(covmat_file, '_block_') + 6  ! index before block number
    IF (istart-7 < 0) STOP 'bad block number'
    iend = INDEX(covmat_file(istart+1:80), '_') ! index afer block number
    ilength = iend-1
    iend = iend + istart

    stemp = covmat_file(istart+1:iend-1)
    READ(stemp, *) blok
    blok = blok + 1
    WRITE (stemp, *) blok

    covmat_file = covmat_file(1:istart) &
         // TRIM(ADJUSTL(stemp)) &
         // covmat_file(iend:80)

  END SUBROUTINE update_covmat_file

END PROGRAM merge_blocks
