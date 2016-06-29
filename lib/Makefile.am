AM_CXXFLAGS=$(OPTS) -I$(top_srcdir)
AM_LDFLAGS=$(LDOPTS)
#if set LDCHOL here instead of in configure.ac compiles "No rule to make target ..."

nobase_include_HEADERS=accphi.h mkg.h readcfg.h mkh.h aos.h mkw.h turbulence.h	\
  mkz.h proj.h laplacian.h pcg.h muv.h draw.h genotf.h fractal.h stfun.h	\
  servo.h hyst.h slaving.h mkdtf.h psd.h kalman.h cn2est.h misc.h locfft.h	\
  zernike.h

EXTRA_DIST=fractal_do.cpp prop_grid.cpp 

AOS_SRC=readcfg.cpp turbulence.cpp mkg.cpp mkh.cpp mkz.cpp mkw.cpp accphi.cpp proj.cpp		\
laplacian.cpp pcg.cpp muv.cpp draw.cpp genotf.cpp fractal.cpp stfun.cpp servo.cpp hyst.cpp	\
slaving.cpp mkdtf.cpp psd.cpp kalman.cpp cn2est.cpp misc.cpp locfft.cpp zernike.cpp

lib_LTLIBRARIES=libaos.la
libaos_la_SOURCES=$(AOS_SRC)
libaos_la_LIBADD=../math/libaomath.la

