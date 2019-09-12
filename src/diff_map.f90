! Simple command line utility to difference fits maps
! Handles correctly the Healpix bad value
!
! 2011-01-11 Reijo Keskitalo

program diff_map

  use fitstools
  use head_fits
  use pix_tools
  use healpix_types
  use extension

  implicit none

  character(len=filenamelen) :: filename_in1, filename_in2, filename_out, argument
  integer(i4b) :: nside, nmaps, ordering, ierr, iarg
  integer(i4b) :: nmaps2, ordering2, nside2
  real(dp), allocatable :: dmap1(:,:), dmap2(:,:)
  real(sp), allocatable :: smap(:,:)
  character(len=80) :: header(1000), header2(1000)
  integer(i8b) :: npixtot, npix
  logical(lgt) :: there, do_double=.false.

  if (narguments() < 2) then
     write (*,'(/,a,/)') ' Usage: diff_map <map1> <map2> [-o <map_out>] [-d]'
     stop
  end if
    
  call getargument(1, filename_in1)
  inquire(file=trim(filename_in1), exist=there)
  if (.not. there) stop 'ERROR: Input map1 does no exist.'
     
  call getargument(2, filename_in2)
  inquire(file=trim(filename_in2), exist=there)
  if (.not. there) stop 'ERROR: Input map2 does no exist.'
     
  iarg = 3
  filename_out = '!diff.fits'
  DO
     IF (iarg > narguments()) EXIT

     CALL getArgument(iarg, argument)
     IF (INDEX(argument, '-d') /= 0) THEN
        do_double = .true.
        write (*,*) ' Producing double precision outputs'
     ELSE IF (INDEX(argument, '-o') /= 0) THEN
        iArg = iArg + 1           
        CALL getArgument(iarg, argument)
        filename_out = TRIM(argument)
     ELSE
        WRITE (*,*) 'Unrecognized argument: '//TRIM(argument)
        STOP
     END IF

     iarg = iarg + 1
  END DO


  header = ''
  header2 = ''

  npixtot = getsize_fits(filename_in1, nmaps=nmaps, ordering=ordering, nside=nside)
  npixtot = getsize_fits(filename_in2, nmaps=nmaps2, ordering=ordering2, nside=nside2)

  if (nside /= nside2) then
     write (*, '(a, i0, a, i0)') ' ERROR: nsides do not match: ', nside, ' != ', nside2
     stop
  end if

  if (nmaps /= nmaps2) then
     write (*, '(a, i0, a, i0)') ' ERROR: nmaps do not match: ', nmaps, ' != ', nmaps2
     stop
  end if

  npix = npixtot
  
  WRITE (*,'(a,i9)')  ' npix     == ', npixtot
  WRITE (*,'(a,i9)')  ' nside    == ', nside
  WRITE (*,'(a,i9)')  ' nmaps    == ', nmaps
  WRITE (*,'(a,i9)')  ' ordering == ', ordering

  ALLOCATE(dmap1(0:npix-1, nmaps), dmap2(0:npix-1, nmaps),  stat=ierr)
  IF (ierr /= 0) STOP 'No room to store the maps!?'

  CALL input_map(filename_in1, dmap1, npix, nmaps, header=header, fmissval=hpx_dbadval)

  CALL input_map(filename_in2, dmap2, npix, nmaps, header=header2, fmissval=hpx_dbadval)

  if (ordering /= ordering2) then
     if (ordering == 1 .and. ordering2 == 2) then
        write (*,'(a)') ' Converting map2 to RING'
        call convert_nest2ring(nside, dmap2)
     else if (ordering == 2 .and. ordering2 == 1) then
        write (*,'(a)') ' Converting map2 to NESTED'
        call convert_ring2nest(nside, dmap2)
     else
        write (*, '(a,i0,a,i0)') 'ERROR: Incompatible orderings: ', ordering, ' != ', ordering2
        stop
     end if
  end if

  where (dmap1 /= hpx_dbadval .and. dmap2 /= hpx_dbadval)
     dmap1 = dmap1 - dmap2
  elsewhere
     dmap1 = hpx_dbadval
  end where
  
  call add_card(header, 'CREATOR', 'diff_map', 'Software creating the file')
  call add_card(header, 'MAP1', trim(filename_in1), 'First map')
  call add_card(header, 'MAP2', trim(filename_in2), 'Second map')
  
  call merge_headers(header, header2)

  if (do_double) then
     call output_map(dmap1, header, filename_out)
  else
     ALLOCATE(smap(0:npix-1, nmaps), stat=ierr)
     IF (ierr /= 0) STOP 'No room to store the map!?'
     smap = real(dmap1, sp)
     call output_map(smap, header, filename_out)
  end if

  write (*,*) 'Stored the difference map into ' // trim(filename_out)

  deallocate(dmap1, dmap2)
  if (.not. do_double) deallocate(smap)

end program diff_map
