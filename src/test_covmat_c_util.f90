program test_covmat_c_util
  implicit none

  external fsize_c

  integer(8) :: sz
  character(len=1024) :: fname='test.dat'

  call fsize_c(fname, len_trim(fname), sz);

  print *,'File size is ',sz,' bytes'
  
end program test_covmat_c_util
