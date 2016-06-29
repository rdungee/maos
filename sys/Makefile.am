AM_CXXFLAGS=$(OPTS) -I$(top_srcdir) 
AM_LDFLAGS=$(LDOPTS)

nobase_include_HEADERS=mem.h process.h scheduler.h sock.h sockio.h thread.h	\
 daemonize.h scheduler_client.h sys.h thread_pool.h misc.h hashlittle.h path.h	\
 common.h bin.h readstr.h
EXTRA_DIST=process_apple.cpp process_linux.cpp process_bsd.cpp process_cygwin.cpp
SYS_SRC=daemonize.cpp misc.cpp process.cpp thread.cpp hashlittle.cpp mem.cpp path.cpp	\
 scheduler_client.cpp sockio.cpp thread_pool.cpp sock.cpp bin.cpp readstr.cpp


lib_LTLIBRARIES=libaosys.la
libaosys_la_SOURCES=$(SYS_SRC)
libaosys_la_LIBADD=-lz -lm -ldl #don't remove
#we put -lz -lm -ldl here so that libtool will always put these flags after the library when linking exes.
