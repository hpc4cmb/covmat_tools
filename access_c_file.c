#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

FILE *matrix_file;
char matrix_filename[255];

void read_c_matrix_(char *filename, double *matrix, int *elements) {
  /* Small routine to read a binary C file */
  int count, i;
  
  /* fix the received filename */
  for (i=0; filename[i] != ' '; i++) ;
  filename[i] = '\0';

  matrix_file = fopen(filename, "r");
  if (!matrix_file) {
    fprintf(stderr, "Unable to open %s!\n", filename);
    return;
  }

  count = fread(matrix, sizeof(double), *elements, matrix_file);
  printf("Succesfully read %i elements from %s\n", count, filename);
  
  fclose(matrix_file);
}

void open_c_matrix_(char *filename, char *mode) {
  /* opens binary matrix file. Use this
     if the matrix is too big to be read all at once
     filename -- pointer to a character string terminated by a blank
                 character
     mode     -- access mode for file (terminated by a blank). Usually
                 "r ", or "w "
   */
  int i;

  /* fix the received filename and mode*/
  for (i=0; filename[i] != ' '; i++) ;
  filename[i] = '\0';

  for (i=0; mode[i] != ' '; i++) ;
  mode[i] = '\0';

  matrix_file = fopen(filename, mode);
  if (!matrix_file) {
    fprintf(stderr, "Unable to open %s!\n", filename);
    return;
  } else printf("Opened %s .\n", filename);

  strcpy(matrix_filename, filename);
}

void close_c_matrix_() {
  if (!fclose(matrix_file)) printf("Matrix file closed succesfully.\n");
  else fprintf(stderr, "Error closing %s.\n", matrix_filename);
}

void read_c_matrix_line_(double *line, int *elements) {
  /* Small routine to read a binary C file */
  int count, i;
  
  count = fread(line, sizeof(double), *elements, matrix_file);
  /*printf("Succesfully read %i elements from %s\n", count, matrix_filename);*/
}


void read_specified_c_matrix_line_(double *line, int *elements, int *lineNum) {
  /* Small routine to read a binary C file */
  int count, i;

  if (fseek(matrix_file,(long)(*elements*sizeof(double) * *lineNum),SEEK_SET)){
    fprintf(stderr, "Unable to find specified location in file.\n");
    return;
  }
  count = fread(line, sizeof(double), *elements, matrix_file);
  /*printf("Succesfully read %i elements from %s\n", count, matrix_filename);*/
}

void write_c_matrix_line_(double *line, int *elements) {
  /* Small routine to write a binary C file */
  int count, i;
  
  count = fwrite(line, sizeof(double), *elements, matrix_file);
  /*printf("Succesfully wrote %i elements to %s\n", count, matrix_filename);*/
}
