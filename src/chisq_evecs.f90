! This code reads madam pixel-pixel inverse covariance matrix and a
! given residual noise map. It then computes the chi squared
! test for the given map using the inverse covariance matrix

PROGRAM chisq_test

  USE healpix_types
  USE fitstools, ONLY   : input_map, getsize_fits
  USE pix_tools, ONLY   : convert_ring2nest
  USE covmat_util, ONLY : tic, toc
  USE extension, ONLY   : getArgument, nArguments

  IMPLICIT NONE

  CHARACTER(len=filenamelen) :: argument, evec_file, eval_file
  character(len=filenamelen) :: mapfile, form
  INTEGER(i4b) :: nstokes, nside, ordering
  integer(i4b) :: isig, iarg, rec_len
  INTEGER(i8b) :: npix, icol
  INTEGER, PARAMETER :: covmat_unit=55
  REAL(dp), POINTER :: evec(:), evals(:), map(:,:), map2(:)
  REAL(dp) :: chisq, alpha, w, evalmax
  REAL(dp) :: dof, limit=1E-30
  LOGICAL :: there, invert=.false.
  INTEGER(i4b) :: ierr, nhit

  IF (nArguments() < 3)  THEN
     WRITE(*,'(/,a,/)') ' Usage: chisq <covmat-evecs> <covmat-evals> ' //&
          '<map> [-inv] [-lim <limit>]'
     STOP
  ELSE
     CALL getArgument(1, argument)
     evec_file = TRIM(argument)
     INQUIRE(file=TRIM(evec_file), exist=there)
     IF (.NOT. there) THEN
        WRITE (*,*) TRIM(evec_file) // ' does not exist!'
        STOP
     END IF
     WRITE(*,*) ' Reading ' // TRIM(evec_file)

     CALL getArgument(2, argument)
     eval_file = TRIM(argument)
     INQUIRE(file=TRIM(eval_file), exist=there)
     IF (.NOT. there) THEN
        WRITE (*,*) TRIM(eval_file) // ' does not exist!'
        STOP
     END IF
     WRITE(*,*) ' Reading ' // TRIM(eval_file)

     CALL getArgument(3, argument)
     mapfile = TRIM(argument)
     INQUIRE(file=mapfile, exist=there)
     IF (.NOT. there) THEN
        WRITE (*,*) TRIM(mapfile) // ' does not exist!'
        STOP
     END IF
     WRITE(*,*) ' Reading ' // TRIM(mapfile)

     iarg = 4
     do
        if (iarg > nArguments()) exit
        CALL getArgument(iarg, argument)

        if (index(argument, '-inv') > 0) then
           invert = .true.
           write (*,*) ' Inverting the read eigenvalues'
        else if (index(argument, '-lim') > 0) then
           iarg = iarg + 1
           CALL getArgument(iarg, argument)
           read(argument, *) limit
           WRITE (*,*) ' omitting eigenmodes having eval/eval0 <',limit
        end if

        iarg = iarg + 1
     end do
  END IF

  ! Read in the map
  CALL tic
  npix = getsize_fits(mapfile, nmaps=nstokes, ordering=ordering, &
       nside=nside)
  ALLOCATE(map(0:npix-1, nstokes), map2(0:npix*nstokes-1), stat=ierr)

  IF (ierr /= 0) STOP 'No room for maps'
  CALL input_map(mapfile, map, npix, nstokes)
  IF (ordering == 1) THEN
     ! Covariance matrix is in NESTED scheme, so must the map be
     WRITE (*,*) ' Note: converting map into NESTED scheme'
     CALL convert_ring2nest(nside, map)
  END IF
  CALL toc('read map')

  ! Remove T monopole
  !monopole = SUM(map(:,1))/npix
  !WRITE (*,'(/,a,ES15.5)') '  Removing temperature monopole == ', monopole
  !map(:,1) = map(:,1) - monopole

  ! write the map into a vector
  do isig = 1, nstokes
     map2((isig-1)*npix:isig*npix-1) = map(:,isig)
  end do

  ! Perform (m^T N^-1 m) by storing only a single eigenvector of the
  ! noise covariance at a time.
  CALL tic
  ALLOCATE(evec(0:npix*nstokes-1), evals(0:npix*nstokes-1), stat=ierr)
  IF (ierr /= 0) STOP 'No room for evecs'

  open(unit=covmat_unit, file=trim(eval_file), status='old', form='formatted')
  read(covmat_unit, *) evals
  close(covmat_unit)

  inquire(iolength=rec_len) evec
  OPEN(unit=covmat_unit, file=TRIM(evec_file), status='old', &
       form='unformatted', access='direct', recl=rec_len)

  open(unit=covmat_unit+1, file=trim(adjustl(mapfile))//'.chisq', &
       status='replace', form='formatted')
  write (*,*) ' Storing the calculation in ',trim(adjustl(mapfile))//'.chisq'
  nhit  = 0; chisq = 0; alpha = 0; w = 0
  evalmax = maxval(evals)
  DO icol = 0, npix*nstokes-1
     alpha = 0
     if (evals(icol)/evalmax > limit) then
        READ(covmat_unit, rec=icol+1) evec
        nhit  = nhit + 1
        alpha = dot_product(map2, evec)
        if (invert) then
           w = alpha**2/evals(icol) 
        else 
           w = alpha**2*evals(icol)
        endif
        chisq = chisq + w
     else
        write (*,'(a,i0)') 'Skipping eigenmode # ', icol
     end if
     write (covmat_unit+1,'(i8,3ES15.5)') icol, w, alpha, chisq
  END DO

  CLOSE(covmat_unit)
  CLOSE(covmat_unit+1)
  CALL toc('compute chi-squared')

  form = '(a,10es15.5)'

  WRITE (*,form) ' Total chi squared ==  ', chisq

  WRITE (*,*)
  WRITE (*,'(i8,a)') nhit, ' included eigenmodes (dof)'
  dof = nhit
  WRITE (*,form) ' chi squared / dof ==  ', chisq/dof
  WRITE (*,*)
  
  WRITE (*,form) ' Expectance == ', REAL(dof,dp)
  WRITE (*,form) '        std == ', (2*dof)**0.5
  WRITE (*,'(/,a,G15.5,a)') &
       ' Deviation from chi squared expected value : ', &
       REAL((chisq-dof),dp)/(2*dof)**0.5, ' sigma'
  WRITE (*,*)
       
END PROGRAM chisq_test
