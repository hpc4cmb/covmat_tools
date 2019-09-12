! Simple command line utility to add fits maps with arbitrary weights.
! Handles correctly the Healpix bad value
!
! 2012-02-14 Reijo Keskitalo
!
! Revisions:
!    2014-06-13 : add support for mixed polarizations. Make map weights optional.

program add_maps

  use fitstools
  use head_fits
  use pix_tools
  use healpix_types
  use extension

  implicit none

  character(len=filenamelen) :: filenames_in(1000), filename_out, argument
  integer(i4b) :: nside, nmaps, ordering, ierr, iarg, imap, nmaps_in, nmaps_read
  integer(i4b) :: nside2, nmaps1, nmaps2, ordering2, j
  real(dp), allocatable :: dmap(:,:), dmap2(:,:)
  real(sp), allocatable :: smap(:,:)
  real(dp) :: factors(1000)
  character(len=80) :: header(1000), header2(1000)
  character(len=8) :: units, key, ttype
  character(len=72) :: comment
  integer(i8b) :: npixtot, npixtot2, npix, i
  logical(lgt) :: there, do_double=.false., t_only=.false., pol=.false.

  if (narguments() < 2) then
     write (*,'(/,a,/)') ' Usage: add_maps <map1> [<factor1>] <map2> [<factor2>] ... ' &
          // '[-o <map_out>] [-d] [-u <units>] [--t_only] [--pol]'
     stop
  end if
    
  iarg = 1
  call getargument(iarg, filenames_in(1))
  inquire(file=trim(filenames_in(1)), exist=there)
  if (.not. there) then
     print *,'Not found: ',trim(filenames_in(1))
     stop 'ERROR: Input map does no exist.'
  end if
  iarg = iarg + 1
     
  call getargument(iarg, argument)
  read(argument, *, iostat=ierr) factors(1)
  if (ierr /= 0 .or. abs(factors(1)) < 1e-30) then
     factors(1) = 1
  else
     iarg = iarg + 1
  end if

  write (*,'(a,g20.10)') ' Map1 = ' // trim(filenames_in(1)) // ', factor1 = ', factors(1)

  nmaps_read = 1
  filename_out = '!sum.fits'
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
     ELSE IF (INDEX(argument, '-t') /= 0) THEN
        t_only = .true.
        pol = .false.
     ELSE IF (INDEX(argument, '-p') /= 0) THEN
        t_only = .false.
        pol = .true.
     ELSE
        nmaps_read = nmaps_read + 1
        filenames_in(nmaps_read) = argument
        inquire(file=trim(filenames_in(nmaps_read)), exist=there)
        if (.not. there) then
           print *,'Not found: ',trim(filenames_in(nmaps_read))
           stop 'ERROR: Input map does no exist.'
        end if
        
        iArg = iArg + 1
        call getargument(iArg, argument)
        read(argument, *, iostat=ierr) factors(nmaps_read)
        if (ierr /= 0 .or. abs(factors(nmaps_read)) < 1e-30 ) then
           iarg = iarg - 1
           factors(nmaps_read) = 1
        end if

        write (*,'(a,i0,a,i0,a,g20.10)') ' Map', nmaps_read, ' = ' &
             // trim(filenames_in(nmaps_read)) // ', factor', nmaps_read, ' = ', factors(nmaps_read)
             
     END IF

     iarg = iarg + 1
  END DO


  header = ''

  npixtot = getsize_fits(filenames_in(1), nmaps=nmaps1, ordering=ordering, nside=nside)
  nmaps = nmaps1
  if (t_only) nmaps = 1
  if (pol) nmaps = 3
  npix = npixtot
  
  WRITE (*,'(a,i9)')  ' npix     == ', npixtot
  WRITE (*,'(a,i9)')  ' nside    == ', nside
  WRITE (*,'(a,i9)')  ' nmaps    == ', nmaps
  WRITE (*,'(a,i9)')  ' ordering == ', ordering

  ALLOCATE(dmap(0:npix-1, nmaps), stat=ierr)
  IF (ierr /= 0) STOP 'No room to store the map!?'

  print *,'Loading ',trim(filenames_in(1))

  dmap = 0
  nmaps_in = min(nmaps, nmaps1)
  CALL input_map(filenames_in(1), dmap, npix, nmaps_in, header=header, &
       fmissval=hpx_dbadval)

  !where (dmap /= hpx_dbadval) dmap = dmap * factors(1)
  do j = 1, nmaps_in
     do i = 0,npix-1
        if (dmap(i,j) /= hpx_dbadval) dmap(i,j) = dmap(i,j) * factors(1)
     end do
  end do

  write(key, '("MAP",i0)') 1
  call add_card(header, key, trim(filenames_in(1)), 'Map')
  write(key, '("FACTOR",i0)') 1
  call add_card(header, key, factors(1), 'Scaling')

  if (nmaps_read > 1) ALLOCATE(dmap2(0:npix-1, nmaps), stat=ierr)
  IF (ierr /= 0) STOP 'No room to store the map2!?'

  do imap = 2, nmaps_read
     print *,'Loading ',trim(filenames_in(imap))

     npixtot2 = getsize_fits(filenames_in(imap), nmaps=nmaps2, ordering=ordering2, nside=nside2)

     if (t_only) nmaps2 = 1

     if (npixtot /= npixtot2 .or. nside /= nside2) then
        write (*,*) 'ERROR: ' // trim(filenames_in(imap)) // ' is incompatible:'
        WRITE (*,'(a,i9)')  ' npix     == ', npixtot2
        WRITE (*,'(a,i9)')  ' nside    == ', nside2
        stop
     end if
     
     dmap2 = 0
     nmaps_in = min(nmaps, nmaps2)
     CALL input_map(filenames_in(imap), dmap2, npix, nmaps_in, header=header2, &
          fmissval=hpx_dbadval)

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
     

     ! PGF95 WHERE segfaults with large arrays. Bad memory management ...
     !where (dmap2 /= hpx_dbadval) dmap2 = dmap2 * factors(imap)
     !where (dmap /= hpx_dbadval .and. dmap2 /= hpx_dbadval)
     !   dmap = dmap + dmap2
     !elsewhere
     !   dmap = hpx_dbadval
     !end where
     do j = 1, nmaps_in
        do i = 0, npix-1
           if (dmap(i, j) == hpx_dbadval .or. dmap2(i, j) == hpx_dbadval) then
              dmap(i, j) = hpx_dbadval
           else
              dmap(i, j) = dmap(i, j) + dmap2(i, j) * factors(imap)
           end if
        end do
     end do

     write(key, '("MAP",i0)') imap
     call add_card(header, key, trim(filenames_in(imap)), 'Map')
     write(key, '("FACTOR",i0)') imap
     call add_card(header, key, factors(imap), 'Scaling')
  end do
  
  if (len_trim(units) > 0) then
     do imap = 1, nmaps
        write(key, '("TUNIT",i1)') imap
        call add_card(header, key, units, 'Map units')
     end do
  end if

  ! Add sensible column names

  do imap = 1, nmaps
     write(key, '("TTYPE",i1)') imap
     write(ttype, '("column",i0.2)') imap
     write(comment, '("label for field ",i1)') imap
     call add_card(header, key, ttype, comment)
  end do

  if (do_double) then
     call output_map(dmap, header, filename_out)
  else
     ALLOCATE(smap(0:npix-1, nmaps), stat=ierr)
     IF (ierr /= 0) STOP 'No room to store the map!?'
     smap = real(dmap, sp)
     call output_map(smap, header, filename_out)
     deallocate(smap)
  end if

  write (*,*) 'Stored the scaled map into ' // trim(filename_out)

  deallocate(dmap)

end program add_maps
