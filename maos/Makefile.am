SUBDIRS = .
AM_CXXFLAGS=$(OPTS) -I $(top_srcdir) #Compiler flags for each target that doesn't specify per target CFLAGS. Supplements CFLAGS
AM_LDFLAGS=$(LDOPTS)
#the libs can not be in AM_LDFLAGS otherwise in mac the binary will not be updated when the libs change (in MAC)#Notice the order. previous dependent on later
nobase_include_HEADERS=common.h maos.h parms.h types.h utils.h setup_powfs.h \
	setup_aper.h setup_recon.h recon.h  sim.h fdpcg.h genseotf.h	\
	mtch.h ahst.h sim_utils.h setup_surf.h save.h mvm_client.h moao.h \
	recon_utils.h pywfs.h
SRCS=maos.cpp setup_parms.cpp utils.cpp setup_powfs.cpp setup_aper.cpp	\
setup_recon.cpp setup_recon_prep.cpp setup_recon_fit.cpp setup_recon_lsr.cpp		\
setup_recon_mvm.cpp sim.cpp recon.cpp fdpcg.cpp wfsgrad.cpp perfevl.cpp wfsints.cpp		\
genseotf.cpp mtch.cpp cachedm.cpp filter.cpp ahst.cpp sim_utils.cpp setup_surf.cpp		\
recon_utils.cpp moao.cpp save.cpp mvm_client.cpp pywfs.cpp

lib_LTLIBRARIES=libmaos.la
libmaos_la_SOURCES=$(SRCS)
libmaos_la_LIBADD=
if USE_CUDA
libmaos_la_LIBADD+=../cuda/libaogpu.la
#if USE_STATIC
#libmaos_la_LIBADD+=../cuda/dlink.o
#endif
endif
libmaos_la_LIBADD+=../lib/libaos.la 
bin_PROGRAMS=../bin/maos 
___bin_maos_SOURCES=main.cpp
___bin_maos_LDADD=libmaos.la  
if USE_CUDA
___bin_maos_LDADD+=$(libmaos_la_LIBADD) $(LIBCUDA)
#if USE_STATIC
#___bin_maos_LDADD+=-lculibos -ldl -lrt
#endif
endif
___bin_maos_LDADD+=$(LDEXE_MAOS) $(LDEXE_FULL) 
