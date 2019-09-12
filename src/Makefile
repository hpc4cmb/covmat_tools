# Temporary fix for buggy Healpix installation in the CMB module
#HEALPIX_INCLUDE = -I/project/projectdirs/planck/software/keskital/$(TARGET)/Healpix/include
#HEALPIX_LINK = -L/project/projectdirs/planck/software/keskital/$(TARGET)/Healpix/lib -lhealpix $(F90OMPFLAGS)

INCLUDE = $(SCALAPACK_INCLUDE) $(HEALPIX_INCLUDE) $(CFITSIO_INCLUDE) \
	$(FFTW3_INCLUDE)
LINK    = $(SCALAPACK_LINK) $(HEALPIX_LINK) $(CFITSIO_LINK) $(FFTW3_LINK)

F90FLAGS += -g $(F90OMPFLAGS)
#F90FLAGS += -O0 -fbounds-check -Wall -pedantic -fbacktrace -fdump-core
#F90FLAGS += -O0 -Mbounds

UTIL    = scalapack_tools.o covmat_util.o fourier_fftw.o partrans.o downgrade.o covmat_c_util.o

ALL_BINARIES = covmat2fits plot_pixel_covariance covmat2block \
	merge_covmat merge_blocks invert_covmat invert_harmonic_covmat \
	covmat2harmonic_mpi harmonic2bias chisq chisq_harmonic \
	downgrade_map evecs2bias merge_interlaced_blocks \
	scale_covmat mask_covmat add_to_diagonal invert_covmat_sparse \
	merge_covmat_mpi merge_covmat_sparse_mpi cc2covmat covmatmul_mpi \
	chisq_mc chisq_evecs downgrade_evecs covmat2block_sparse \
	generate_noisemap evecs2covmat evecs2sqrtcovmat \
	generate_whitenoisemap merge_blocks_mpi project_bad_eigenmode \
	remove_temperature plot_pixel_covariance_sparse covmat2block_mpi \
	mapstats downgrade_covmat_diagonal evecs2covmat_mpi invert_covmat_mpi \
	downgrade_evecs_mpi nobs2rcond evecs2bias_mpi covmat2commander \
	extract_diagonal generate_noisemap_mpi scale_map diff_map add_maps map2mask add_covmat_mpi

all : $(UTIL) $(ALL_BINARIES)

install : all
	mkdir -p $(PREFIX)/covmat_tools/bin
	cp $(ALL_BINARIES) $(PREFIX)/covmat_tools/bin

% : $(UTIL) %.o
	$(MPF90) $(F90FLAGS) -o $* $*.o $(UTIL) $(LINK) $(LDFLAGS)

%.o : %.F90
	$(MPF90) $(F90FLAGS) $(INCLUDE) -c $<

%.o : %.f90
	$(MPF90) $(F90FLAGS) $(INCLUDE) -c $<

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

clean : 
	$(RM) -rf *~ *.o *.mod *.dSYM $(ALL_BINARIES)
