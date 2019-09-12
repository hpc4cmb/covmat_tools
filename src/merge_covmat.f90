! This program is used to combine detector covariance matrices
! into a single matrix (as inverses they are additive).
!
! February 18th 2008 Reijo Keskitalo
! 
! Revisions :
! July 25th 2008 -- can now merge matrices with different nstokes parameters

PROGRAM merge_covmat

  USE healpix_types
  USE extension, ONLY   : nArguments, getArgument
  USE covmat_util, ONLY : tic, toc

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, covmat_file(100)
  CHARACTER(len=filenamelen) :: outfilename
  INTEGER(i4b) :: nstokes(100), nside, ncovmat, nstokesmax, iArgument, &
       iCovmat, rec_len, iErr
  INTEGER(i8b) :: npix, icol
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: covmatIn(:), covmatOut(:)
  LOGICAL :: there, interlaced=.false.

  IF (nArguments() < 5)  THEN
     write (*,*) 'Usage: merge_covmat <nside> <nstokes1> <covmat1> ' // &
          '[<nstokes2> <covmat2> [<nstokes3> <covmat3> ...]] ' // &
          '[--interlaced] [-o <covmat_out>]'
     stop
  ELSE
     CALL getArgument(1, argument)
     READ (argument, *) nSide
     WRITE (*,'(a,i4)') 'Nside   == ', nSide

     covmat_file = ''
     nStokes     = 0
     nCovmat     = 0
     outfilename = 'merged_covmat.dat'
     iArgument   = 2
     do
        CALL getArgument(iArgument, argument)

        if (index(argument,'-o') /= 0) then
           CALL getArgument(iArgument+1, argument)
           outfilename = TRIM(ADJUSTL(argument))
           iArgument = iArgument + 2
        else if (index(argument,'-inter') /= 0) then
           interlaced = .true.
           WRITE (*,'(a)') 'Treating matrices as INTERLACED'
           iArgument = iArgument + 1
        else
           nCovmat = nCovmat + 1
           READ (argument, *) nStokes(nCovmat)
           WRITE (*,'(a,i0,a,i4)') 'Nstokes', nCovmat, '  == ',nStokes(nCovmat)
           
           CALL getArgument(iArgument+1, argument)
           covmat_file(nCovmat) = TRIM(argument)
           INQUIRE(file=TRIM(covmat_file(nCovmat)), exist=there)
           IF (.NOT. there) STOP 'covmat_file does not exist!'
           WRITE (*,'(a,i0,a)') &
                'covmat', nCovmat, '  == ' // trim(covmat_file(nCovmat))

           iArgument = iArgument + 2
        end if
        if (iArgument > nArguments()) exit
     end do

     WRITE (*,*) ' Writing to ' // TRIM(outfilename)
  END IF

  nPix = 12*nSide**2
  nStokesMax = maxval(nStokes)

  ALLOCATE(covmatIn(0:npix*nstokesmax-1), covmatOut(0:npix*nstokesmax-1), &
       stat=iErr)
  IF (iErr /= 0) STOP 'no room for covmat'

  do iCovmat = 1, ncovmat
     inquire(iolength=rec_len) covmatIn(0:nPix*nStokes(iCovmat)-1)
     OPEN(unit=covmat_unit+iCovmat, file=TRIM(covmat_file(iCovmat)), &
          status='old', form='unformatted', access='direct', &
          recl=rec_len)
  end do

  inquire(iolength=rec_len) covmatIn(0:nPix*nStokesMax-1)
  OPEN(unit=covmat_unit, file=TRIM(ADJUSTL(outfilename)), status='replace', &
       form='unformatted', access='direct', recl=rec_len)

  CALL tic
  DO icol = 0, nPix*nStokesMax - 1
     covmatOut = 0.0
     do iCovmat = 1, nCovmat
        if (.not. interlaced .and. icol >= nPix*nStokes(iCovmat)) cycle
        if (interlaced .and. nStokes(iCovmat) < nStokesMax) then
           ! merging unpolarized matrix to interlaced, polarized matrix
           if (modulo(icol, int(nStokesMax, i8b)) /= 0) cycle
           
           READ (covmat_unit+iCovmat, rec=icol/nStokesMax+1) &
                covmatIn(0:nPix*nStokes(iCovmat)-1)
           
           ! replace sentinel value -1 by zero for merging
           if ( covmatIn(icol) < 0 ) covmatIn(icol) = 0
                
           covmatOut(0:nPix*nStokesMax-1:nStokesMax) = &
                covmatOut(0:nPix*nStokesMax-1:nStokesMax) &
                + covmatIn(0:nPix*nStokes(iCovmat)-1)
        else
           ! all other cases
           READ (covmat_unit+iCovmat, rec=icol+1) covmatIn(0:nPix*nStokes(iCovmat)-1)

           ! replace sentinel value -1 by zero for merging
           if ( covmatIn(icol) < 0 ) covmatIn(icol) = 0
                
           covmatOut(0:nPix*nStokes(iCovmat)-1) = covmatOut(0:nPix*nStokes(iCovmat)-1) &
                + covmatIn(0:nPix*nStokes(iCovmat)-1)
        end if
     end do

     ! re-introduce the sentinel value
     if ( covmatOut(icol) == 0 ) covmatOut(icol) = -1

     WRITE (covmat_unit, rec=icol+1) covmatOut
  END DO
  CALL toc('merge covmat')

  do icovmat = 0, ncovmat
     CLOSE(unit=covmat_unit+icovmat)
  end do

  DEALLOCATE(covmatIn, covmatOut)

END PROGRAM merge_covmat
