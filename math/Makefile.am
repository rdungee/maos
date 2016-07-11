AM_CXXFLAGS=$(OPTS) -I$(top_srcdir)
AM_LDFLAGS=$(LDOPTS)

nobase_include_HEADERS= type.h numtype.h mathmisc.h random.h mathdef.h cell.h  suitesparse.h blas.h \
 matbin.h spbin.h zfarr.h chol.h  mat.h sp.h fft.h matmath.h loc.h locbin.h defs.h

EXTRA_DIST=suitesparse.cpp matmmap.cpp
EXTRA_DIST+=mat.cpp matmath.cpp blas.cpp matbin.cpp fft.cpp  matcomp.cpp amoeba.cpp bessik.cpp sp.cpp spmm.cpp spbin.cpp 

AOMATH_SRC= dmath.cpp cmath.cpp smath.cpp zmath.cpp lmath.cpp mathmisc.cpp random.cpp zfarr.cpp chol.cpp loc.cpp locbin.cpp 

lib_LTLIBRARIES=libaomath.la
libaomath_la_SOURCES=$(AOMATH_SRC)
libaomath_la_LIBADD=../sys/libaosys.la $(LDCHOL)
# $(LDFFT) $(LAPACK) don't put fft or lapack here for flexibility when compiling for mex routines
