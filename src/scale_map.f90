! Simple command line utility to scale a fits map.
! Handles correctly the Healpix bad value
!
! 2011-01-11 Reijo Keskitalo

program scale_map

  use fitstools
  use head_fits
  use pix_tools
  use healpix_types
  use extension

  implicit none

  character(len=filenamelen) :: filename_in, filename_out, argument
  integer(i4b) :: nside, nmaps, ordering, ierr, iarg, imap
  real(dp), allocatable :: dmap(:, :)
  real(sp), allocatable :: smap(:, :)
  real(dp) :: factor
  character(len=80) :: header(1000)
  character(len=8) :: units, key
  integer(i8b) :: npixtot, npix
  logical(lgt) :: there, do_double=.false.

  if (narguments() < 2) then
     write (*,'(/,a,/)') ' Usage: scale_map <map> <factor> [-o <map_out>] [-d] [-u <units>]'
     stop
  end if
    
  call getargument(1, filename_in)
  inquire(file=trim(filename_in), exist=there)
  if (.not. there) stop 'ERROR: Input map does no exist.'
     
  call getargument(2, argument)
  read(argument, *, iostat=ierr) factor
  if (ierr /= 0) then
     write (*,*) 'ERROR: Failed to parse ', trim(argument), ' for scaling.'
     stop
  end if

  iarg = 3
  filename_out = '!' // trim(filename_in)
  units = ''
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
     ELSE IF (INDEX(argument, '-u') /= 0) THEN
        iArg = iArg + 1           
        CALL getArgument(iarg, argument)
        units = TRIM(argument)
     ELSE
        WRITE (*,*) 'Unrecognized argument: '//TRIM(argument)
        STOP
     END IF

     iarg = iarg + 1
  END DO


  header = ''

  npixtot = getsize_fits(filename_in, nmaps=nmaps, ordering=ordering, nside=nside)
  npix = npixtot
  
  WRITE (*,'(a,i9)')  ' npix     == ', npixtot
  WRITE (*,'(a,i9)')  ' nside    == ', nside
  WRITE (*,'(a,i9)')  ' nmaps    == ', nmaps
  WRITE (*,'(a,i9)')  ' ordering == ', ordering

  ALLOCATE(dmap(0:npix-1, nmaps), stat=ierr)
  IF (ierr /= 0) STOP 'No room to store the map!?'

  CALL input_map(filename_in, dmap, npix, nmaps, header=header, &
       fmissval=hpx_dbadval)

  where (dmap /= hpx_dbadval) dmap = dmap * factor
  
  call add_card(header, 'factor', factor, 'Scaled by scale_map')

  if (len_trim(units) > 0) then
     do imap = 1, nmaps
        write(key, '("TUNIT",i1)') imap
        call add_card(header, key, units, 'Map units')
     end do
  end if


  if (do_double) then
     call output_map(dmap, header, filename_out)
  else
     ALLOCATE(smap(0:npix-1, nmaps), stat=ierr)
     IF (ierr /= 0) STOP 'No room to store the map!?'
     smap = real(dmap, sp)
     call output_map(smap, header, filename_out)
  end if

  write (*,*) 'Stored the scaled map into ' // trim(filename_out)

  deallocate(dmap)
  if (.not. do_double) deallocate(smap)

end program scale_map
