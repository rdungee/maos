AM_CXXFLAGS=$(MEXOPT) -I$(SRC_DIR)/lib -I$(SRC_DIR)/math 
AM_LDFLAGS=$(LDMEX)
#LDFLAGS goes before LDADD. Must put maos libraries before matlab libraries.
#The following should be in LDADD instead of AM_LDFLAGS to recompile exes when they are updated
LDADD=../lib/libaos.la $(LIBMEX) 
EXTRA_DIST=io.cpp io.h interface.h Makefile.mex compile.sh aolib.h aolibmex.h
EXEEXT=.$(MEXEXT) #Automatic add suffix

.NOTPARALLEL:io.cpp

bin_PROGRAMS=read write aolib maos 

#.cpp.o:
#	$(MEX) $(AM_CXXFLAGS) -c $< 

#Two special ones that does not depend on aos library
#read$(EXEEXT):read.cpp io.cpp
#	$(MEX) $(SRC_DIR)/mex/read.cpp $(SRC_DIR)/mex/io.cpp -lz
#write$(EXEEXT):write.cpp io.cpp
#	$(MEX) $(SRC_DIR)/mex/write.cpp $(SRC_DIR)/mex/io.cpp -lz

read_SOURCES=read.cpp io.cpp
read_CXXFLAGS=$(MEXOPT)
read_LDADD=-lmx -lmex -lz

write_SOURCES=write.cpp io.cpp
write_CXXFLAGS=$(MEXOPT)
write_LDADD=-lmx -lmex -lz

maos_LDADD=../maos/libmaos.la 
if USE_CUDA
maos_LDADD+=$(LIBCUDA) #must be before LIBMEX, otherwise linking error in debug mode
endif
maos_LDADD+=-lcholmod $(LIBMEX) 


BUILT_SOURCES=maos2mex.h aolib.cpp aolibmex.cpp
maos2mex.h:../maos/parms.h ../maos/types.h ../lib/aos.h ../scripts/maos2mex.py
	$(SRC_DIR)/scripts/maos2mex.py $(SRC_DIR) $(SRC_DIR)/mex/maos2mex.h
aolib.cpp:../maos/parms.h ../maos/types.h ../scripts/aolib.py aolib.h
	$(SRC_DIR)/scripts/aolib.py $(SRC_DIR) $(SRC_DIR)/mex/aolib.cpp
aolib_SOURCES=aolibmex.cpp aolib.cpp
maos_SOURCES=maos.cpp
