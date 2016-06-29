SUBDIRS = . 
AM_CXXFLAGS=$(OPTS) -I $(top_srcdir)
AM_LDFLAGS=$(LDOPTS)

bin_PROGRAMS=../bin/skyc 
___bin_skyc_SOURCES=parms.h types.h setup_parms.cpp  skyc.cpp  skyc.h  utils.cpp  utils.h skysim.cpp skysim.h setup_powfs.cpp setup_powfs.h setup_aster.cpp setup_aster.h photon.cpp photon.h skysim_utils.cpp skysim_utils.h mtch.cpp mtch.h  nafocus.cpp nafocus.h genstars.cpp genstars.h genpistat.cpp genpistat.h setup_star.cpp setup_star.h
___bin_skyc_LDADD=../lib/libaos.la $(LDEXE_FULL)
