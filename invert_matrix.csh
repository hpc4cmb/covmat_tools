#!/bin/csh -f

if ($#argv < 2) then
  echo "Usage: $0 <matrix> <ntasks>"
else
  set matrixfile = $1
  set ntasks = $2
  echo "inverting $matrixfile using $ntasks cores"

  # compute size of the matrix
  set bytes = `ls -l $matrixfile | awk '{print $5}'`
  set npix  = `echo "${bytes}/8.0" | bc -l`
  set npix  = `printf "%15.0f" $npix`
  set nrows = `echo "sqrt(${bytes}/8.0)" | bc -l`
  set nrows = `printf "%15.0f" $nrows`
  echo "Matrix has $npix pixels and $nrows rows"

  # compute number of cores per row
  set ncores_per_row = `echo "sqrt(${ntasks})" | bc -l`
  set ncores_per_row = `printf "%15.0f" $ncores_per_row`
  echo "Distributing for $ncores_per_row processor rows"

  # block distribute the matrix
  mkdir -p files
  cd files
  if (-e $matrixfile) then
  else
    ln -s ../$matrixfile
  endif
    
  cmb_blockdist $matrixfile $nrows $nrows $ncores_per_row $ncores_per_row 32 1

  # create a pbs script for inverting the matrix
  cd ..
  set jobfile = invert_${matrixfile}.pbs
  rm -f $jobfile
  touch $jobfile
  echo "#PBS -q debug" >> $jobfile
  echo "#PBS -l mppwidth=$ntasks" >> $jobfile
  echo "#PBS -l walltime=00:30:00" >> $jobfile
  echo "#PBS -j eo" >> $jobfile
  echo "#PBS -V" >> $jobfile
  echo "" >> $jobfile
  echo 'cd $PBS_O_WORKDIR' >> $jobfile
  echo "" >> $jobfile
  echo "aprun -n $ntasks mad_ppeigen -dist $nrows -evecs -block 32 $matrixfile" >> $jobfile
  echo "aprun -n $ntasks mad_ppinvert -dist $nrows -reghigh 1e-6 -block 32 $matrixfile" >> $jobfile
#  echo "aprun -n $ntasks mad_ppinvert -dist $nrows -block 32 $matrixfile" >> $jobfile

  echo "#" | tee -a $jobfile
  echo "# Please submit the job using" | tee -a $jobfile
  echo "#    qsub $jobfile" | tee -a $jobfile
  echo "# Then recompose the matrix using" | tee -a $jobfile
  echo "#   cmb_blockdist ${matrixfile}-evecs $nrows $nrows $ncores_per_row $ncores_per_row 32 -1" | tee -a $jobfile
  echo "#   cmb_blockdist ${matrixfile}.inverse $nrows $nrows $ncores_per_row $ncores_per_row 32 -1" | tee -a $jobfile
  echo "# in the 'files' subdirectory" | tee -a $jobfile
  echo "#" | tee -a $jobfile
endif
