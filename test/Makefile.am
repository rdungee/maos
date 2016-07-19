AM_CXXFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS) 
#LDADD=../lib/libaos.la $(LDEXE_FULL)
AM_DEFAULT_SOURCE_EXT=.cpp
#bin_PROGRAMS=test_mat test_plot test_draw test_rand test_sparse test_cmat2		\
#		test_fft test_proj test_loc test_chol test_accphi convert_psfout	\
#		test_fd test_atm test_fractal test_svdlas2 test_fits 
#if USE_CUDA
#bin_PROGRAMS+=test_gpu
#test_gpu_LDADD=../cuda/libaogpu.la $(LDADD)  $(LIBCUDA)
#endif
LDADD=../sys/libaosys.la $(LDEXE_FULL)
nobase_include_HEADERS=catch.h
bin_PROGRAMS=test
test_SOURCES=test.cpp test_array.cpp
