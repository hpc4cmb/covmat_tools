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

  CHARACTER(len=80)      :: argument, covmat_file, first_file, outfilename
  INTEGER(dp)            :: npix, nstokes, nside
  INTEGER, PARAMETER     :: covmat_unit=55
  REAL(dp), POINTER      :: covmat_in(:,:), covmat_out(:,:)
  LOGICAL                :: there
  INTEGER(i4b)           :: ierr, icol_in, icol_out, n, rec_len

  IF (nArguments() /= 3)  THEN
     STOP 'Usage: merge_blocks <first_block> <nside> <nstokes>'
  ELSE
     CALL getArgument(1, argument)
     first_file = TRIM(argument)
     INQUIRE(file=TRIM(first_file), exist=there)
     IF (.NOT. there) STOP 'covmat_file does not exist!'
     WRITE (*,*) ' Reading ' // TRIM(first_file)

     CALL getArgument(2, argument)
     READ (argument, *) nside
     WRITE (*,'(a,i4)') 'Nside   == ', nside

     CALL getArgument(3, argument)
     READ (argument, *) nstokes
     WRITE (*,'(a,i4)') 'Nstokes == ', nstokes
  END IF

  CALL get_outfilename(first_file, outfilename)
  WRITE (*,*) 'Writing to ' // TRIM(outfilename)

  npix = 12*nside**2
  n    = 1
  if (nstokes == 3) n = 9

  ALLOCATE(covmat_in(n, 0:npix-1), covmat_out(0:npix*nstokes-1, nstokes), &
       stat=ierr) ! need only one column
  IF (ierr /= 0) STOP 'no room for covmat'

  CALL tic
  inquire(iolength=rec_len) covmat_out
  OPEN(unit=covmat_unit+1, file=TRIM(ADJUSTL(outfilename)), status='replace',&
       form='unformatted', access='direct', recl=rec_len)

  icol_in     = 1
  icol_out    = 1
  covmat_file = first_file
  loop_files : DO 
     ! cycle through available files reading the same column from
     ! all of the files and writing it into a single file
     INQUIRE(file=TRIM(covmat_file), exist=there)
     IF (there) then
        inquire(iolength=rec_len) covmat_in
        OPEN(unit=covmat_unit, file=TRIM(covmat_file), status='old',&
             form='unformatted', access='direct', recl=rec_len)
        READ (covmat_unit, rec=icol_in) covmat_in
        close(unit=covmat_unit)

        if (nstokes == 1) then
           covmat_out(:,1) = covmat_in(1,:)
        else
           covmat_out(0:npix*nstokes-1:nstokes, 1) = covmat_in(1, :)
           covmat_out(1:npix*nstokes-1:nstokes, 1) = covmat_in(2, :)
           covmat_out(2:npix*nstokes-1:nstokes, 1) = covmat_in(3, :)

           covmat_out(0:npix*nstokes-1:nstokes, 2) = covmat_in(4, :)
           covmat_out(1:npix*nstokes-1:nstokes, 2) = covmat_in(5, :)
           covmat_out(2:npix*nstokes-1:nstokes, 2) = covmat_in(6, :)

           covmat_out(0:npix*nstokes-1:nstokes, 3) = covmat_in(7, :)
           covmat_out(1:npix*nstokes-1:nstokes, 3) = covmat_in(8, :)
           covmat_out(2:npix*nstokes-1:nstokes, 3) = covmat_in(9, :)
        end if

        WRITE (covmat_unit+1, rec=icol_out) covmat_out
        icol_out = icol_out + 1
        if (icol_out > npix) exit loop_files

        CALL update_covmat_file(covmat_file)
     else
        covmat_file = first_file
        icol_in     = icol_in + 1
     end IF
  END DO loop_files

  CLOSE(unit=covmat_unit+1)

  CALL toc('merge blocks')



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
