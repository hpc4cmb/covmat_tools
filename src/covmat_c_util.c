#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void fsize_c_(char *fname, int *fnamelen, long *nbyte) {
  char *my_fname;
  struct stat st;
  int err;

  my_fname = (char*) malloc((*fnamelen+1) * sizeof(char));
  if (!my_fname) {
    *nbyte = -2;
    return;
  }

  strncpy(my_fname, fname, *fnamelen);
  my_fname[*fnamelen] = '\0';

  err = stat(my_fname, &st);

  if (!err) 
    *nbyte = st.st_size;  
  else
    *nbyte = -1;

  free(my_fname);
}
