program use_read_c_file
  real*8, allocatable :: matrix(:,:)
  integer, parameter :: nside=8, nstokes=3
  integer :: npix, N, ierr
  character(len=255) :: filename

  filename = 'test.dat'

  npix=12*nside**2
  N=(nstokes*npix)**2

  allocate(matrix(0:npix*nstokes-1, 0:npix*nstokes-1), stat=ierr)
  if (ierr/=0) stop 'no room for matrix'

  call read_c_matrix(filename, matrix, N)

  write (*,*) matrix(0,0:10)

  call open_c_matrix(filename)
  call read_c_matrix_line(matrix(:,0), npix*nstokes)
  call close_c_matrix()

  write (*,*) matrix(0:10, 0)

end program use_read_c_file
