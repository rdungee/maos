ACLOCAL_AMFLAGS=-I m4 --install
SUBDIRS=
SUBDIRS+=sys math lib
if USE_CUDA
SUBDIRS+=cuda
endif
SUBDIRS+=maos skyc
if USE_MEX
SUBDIRS+=mex
endif
SUBDIRS+=tools test 
include $(top_srcdir)/aminclude.am

#dist-hook:
#	rm -rf $(distdir)/configure.ac 

all-local:$(srcdir)/math/tmp.cpp
$(srcdir)/math/tmp.cpp:$(srcdir)/math/mat.cpp $(srcdir)/math/cell.cpp $(srcdir)/math/matbin.cpp \
	$(srcdir)/math/spbin.cpp $(srcdir)/math/sp.cpp
	@echo "#define USE_COMPLEX" > $(srcdir)/math/tmp.cpp
	@cat $(srcdir)/math/mat.cpp $(srcdir)/math/cell.cpp $(srcdir)/math/matbin.cpp $(srcdir)/math/spbin.cpp $(srcdir)/math/sp.cpp >> $(srcdir)/math/tmp.cpp

EXTRA_DIST = doxygen.cfg docs version.h scripts config
