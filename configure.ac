AC_PREREQ(2.57)
AC_INIT(maos, 2.3.0, lianqiw@gmail.com)
AM_INIT_AUTOMAKE([1.9.6 gnu -Wall no-define subdir-objects])
#silent building
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])]) 
#The URL where misc supportive libaries are located.
BASEURLS=
BASEURLS+=" https://github.com/lianqiw/files/raw/master/"
BASEURLS+=" https://sourceforge.net/projects/lianqiw/files/external"
BASEURLS+=" https://bitbucket.org/lianqiw/files/downloads/"
#If user does not supply any CFLAGS, set something so that AC_PROG_CC won't put a -g -O2 to CFLAGS
if test -z "$CFLAGS" ;then
    CFLAGS=" "
fi
if test -z "$LIBS" ;then
    LIBS=" "
fi
AC_CONFIG_MACRO_DIR([m4])
#After enable HEADER, the package name will not appear in CFLAGS. Appear in config.h instead.
AC_CONFIG_HEADER([config.h])
#defines host host_os host_cpu host_vendor etc
AC_CANONICAL_HOST
AC_LANG([C])
AC_TYPE_SIZE_T #it sets CC if empty
AC_PROG_CC([$CC])
AM_PROG_CC_C_O
#AC_PROG_CC_STDC #enable the current c standard, which is c99 circa 2014
AC_PROG_CC_C99  #enable C99
#AC_PROG_CXX([$CXX])
AC_OPENMP
#The following must be after AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_DISABLE_SHARED
AC_ENABLE_STATIC
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL
#Doxygen Features
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN($PACKAGE_NAME, doxygen.cfg)
if test x$DX_FLAG_doc = x0 ;then
    AM_CONDITIONAL(DO_DOXYGEN, false)
else
    AM_CONDITIONAL(DO_DOXYGEN, true)
fi

AC_C_BIGENDIAN([BIGENDIAN=1], [BIGENDIAN=0], AC_MSG_ERROR("Unknown Endianness"))
AC_DEFINE_UNQUOTED(BIGENDIAN, [$BIGENDIAN], ["Endianness"])

#top_srcdir and abs_top_srcdir is empty. no idea why. but srcdir sometimes is the full path.
SRC_DIR=$(cd "$srcdir" && pwd)
BUILD_DIR=$(pwd)
#AC_DEFINE_UNQUOTED means expand the $ variable
AC_DEFINE_UNQUOTED(SRCDIR,  ["$SRC_DIR"], ["Define SRC folder"])
AC_DEFINE_UNQUOTED(BUILDDIR,["$BUILD_DIR"],["Define BUILD folder"])

#find the real flavor of CC
if $CC -v 2>&1 |grep clang ;then
    CN=CLANG
elif $CC -v 2>&1 |grep icc ;then
    CN=ICC
else
    CN=GCC
fi
#canonicalize the host_cpu name
case "$host_cpu" in
    *86)
	host_cpu=x86
	;;
    *64)
	host_cpu=x86_64
	;;
esac
#AUTO detect libraries in the system
#Try to discover CUDA install in the system
if ! which which >/dev/null 2>&1;then
    die "Please install 'which'"
fi
if test -n "`which nvcc 2>/dev/null`" ;then
    NVCC="$(which nvcc)"
    if test -L "$NVCC" ;then
	NVCC="$(readlink "$NVCC")"
    fi
else
    for cuda in /opt/cuda/default /opt/cuda /usr/local/cuda ;do
	if test -x "${cuda}/bin/nvcc" ;then
	    NVCC="${cuda}/bin/nvcc"
	    break;
	fi
    done
    unset nvcc
fi
if test -n "$NVCC" ;then
    CUDADIR_AUTO="$(dirname "$(dirname "${NVCC}")")"
    unset NVCC
fi

#Compile mex routines if MATLAB is detected and single is not enabled
MEXDIR_AUTO=
if test -n "`which mex 2>/dev/null`" ;then
    MEX="$(which mex)"
    if test -L "$MEX" ;then
	MEX="$(readlink "$MEX")"
    fi
    MEXDIR_AUTO="$(dirname "$MEX")"
else
    for i in /Applications/MATLAB*.app ~/Applications/MATLAB*.app /opt/MATLAB/* /usr/local/share/MATLAB/* ;do
	if test -x "$i/bin/mex" ;then
       	    MEXDIR_AUTO="${i}/bin"
	    echo "Matlab found at ${i}"
	    break;
	fi
    done
fi

#Common configure options.
AC_ARG_ENABLE(all-static, AS_HELP_STRING([--diable-all-static], [Compile statically]), [use_static="$enableval"],[use_static="no"])
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[Turn on debug]), [debug="$enableval"], [debug="no"])
AC_ARG_ENABLE(oprof, AS_HELP_STRING([--enable-oprof],[Enable oprof for code profiling.]),[use_oprof="$enableval"],[use_oprf="no"])
AC_ARG_ENABLE(cuda-double, AS_HELP_STRING([--enable-cuda-double],[Enable cuda double precision]), [cuda_double="$enableval"],[cuda_double="no"])
AC_ARG_ENABLE(scheduler,AS_HELP_STRING([--disable-scheduler],["Disable built in scheduler"]), [use_scheduler="$enableval"], [use_scheduler="yes"])
AC_ARG_ENABLE(websockets,AS_HELP_STRING([--disable-websockets],[Disable websockets (default: enabled)]), [use_websockets="$enableval"], [use_websockets="yes"])
AC_ARG_ENABLE(mkl, AS_HELP_STRING([--enable-mkl], ["Prefer to use MKL"]), [use_mkl="$enableval"], [use_mkl="default"])
AC_ARG_WITH(cuda, AS_HELP_STRING([--with-cuda[=DIR]],[With cuda support]), [CUDADIR=${withval}], [CUDADIR="$CUDADIR_AUTO"])
AC_ARG_WITH(ccbin, AS_HELP_STRING([--with-ccbin=CC,],[Specify CCBIN for nvcc]), [CCBIN=${withval}], [CCBIN=${CC}])
AC_ARG_WITH(gtk,  AS_HELP_STRING([--with-gtk[=DIR]],[Enable gtk support]), [GTKDIR=${withval}], [GTKDIR=])
AC_ARG_ENABLE(gtk-3, AS_HELP_STRING([--enable-gtk-3],[Using GTK+-3]), [use_gtk_3="$enableval"], [use_gtk_3="default"])
AC_ARG_WITH(matlab, AS_HELP_STRING([--with-matlab[=DIR]],[Enable compiling mex routines for matlab]), [MEXDIR=${withval}], [MEXDIR=${MEXDIR_AUTO}])
AC_ARG_ENABLE(release, AS_HELP_STRING([--enable-release], [Compile for broad compatibility]), [use_release="$enableval"],[use_release="no"])

if test -z "$MEXDIR" -o "$MEXDIR" = "no" ;then
    with_matlab=no
else
    with_matlab=yes
fi

#Specifying both -O3 and -g in mac produces internal error in gcc 4.2
#-all-static will produce a static exe if all libraries are available in static linking mode.
OPTScommon=
if test "${USER}" = "lianqiw" ;then
    OPTScommon+=" -Wshadow -Wall -Wextra -Werror -Wno-unknown-pragmas -Wno-unused-function -fno-common -Wno-missing-field-initializers -Wno-format-security -Wno-vla" 
fi
LDEXE=
OPTS=
#Set optimize flags for different compilers
#-xO Enables SSE3, SSE2 and SSE instruction sets optimizations for non-Intel CPUs 
COMMON_ICC="-wd193 -wd1418 -wd2259 -wd1372 -wd981 -wd1292 -wd10237 -wd858 -wd1478 -wd1786 -wd3180 -diag-disable remark -vec-report=3"
COMMON_GCC=" -Wno-deprecated-declarations -fms-extensions -static-libgcc"
COMMON_CLANG="-Wno-deprecated-declarations -Qunused-arguments "
COMMON_ALL=""

OPTIM_ALL="-O3 "
DEBUG_ALL=" -g -O0 -DDEBUG_TRACE"

DEBUG_GCC="-ggdb"

#inline definition:

#in GCC, extern inline will use the inline and not emit function definition.
#in GCC, inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header.
#in C99, extern inline will cause use the inline, but also emit function deminition, causing multiple definitions if put in header. This is equivalent to GCC inline
#in C99, inline will use the inline and not emit function definition. This is equivalent to GCC extern inline 
#GCC mode can be forced on using __attribute__((gnu_inline))
#always_inline causes inline to be done in debug mode.
case "$CN" in
    ICC)
	CC_VERSION=$($CC -dumpversion)
	if test -n "$CC_VERSION" ;then
	    major=${CC_VERSION%%.*}
	    minor=${CC_VERSION#*.}
	    minor=${minor%%.*}
	    if test $major -lt 11 ;then
		AC_MSG_ERROR([OpenMP is not supported by ICC $CC_VERSION])
		unset OPENMP_CFLAGS
	    fi
	fi
	;;
    GCC)
	#Test GCC version
	CC_VERSION=$($CC -dumpversion)
	if test -n "$CC_VERSION" ;then
	    major=${CC_VERSION%%.*}
	    minor=${CC_VERSION#*.}
	    if test $minor = $CC_VERSION ;then
		minor=0
	    else
		minor=${minor%%.*}
	    fi
	    echo major=$major minor=$minor
	    if test "$major" -lt 4 ;then
		AC_MSG_ERROR([At least version 4 is required when compiling with gcc.])
	    fi
	    OPTIM_GCC+=" -ftree-vectorize " #introduced in GCC 4
	    if test "$major" -gt 4 -o "$minor" -ge 2 ;then
		#echo "Enabling -mtune=native"
		OPTIM_GCC+=" -mtune=native" #introduced in GCC 4.2
	    fi
	    if test "$major" -ge 7 ;then #GCC 7 introduces more checking.
		COMMON_GCC+=" -Wno-format-truncation -Wno-misleading-indentation -Wimplicit-fallthrough=2"
	    fi
	    #GCC4.3 changed meaning of inline to C99
	    #if test $major -eq 4 -a $minor -lt 5;then
	    #gnu_inline="yes"
	    #fi
	    #GCC introduced openmp task in 4.4
	    if test $major -le 4 -a $minor -le 4 ;then
		AC_MSG_ERROR([OpenMP is not supported by GCC $CC_VERSION])
		unset OPENMP_CFLAGS
	    fi
	fi
	;;
    CLANG)
	if test "$system" = "apple" ;then
	    CFLAGS+=" -stdlib=libc++"
	fi
	;;
    *)
	;;	
esac
if test "x$debug" = "xyes" ;then
    AM_CONDITIONAL(DEBUG,true)
    AC_DEFINE(DEBUG,[1],["Debug the code"])
    ccprefix=DEBUG
else
    AM_CONDITIONAL(DEBUG,false)
    AC_DEFINE(DEBUG,[0],["Do not debug"])
    ccprefix=OPTIM
fi
optname1=COMMON_ALL
optname2=COMMON_${CN}
optname3=${ccprefix}_ALL
optname4=${ccprefix}_${CN}
echo OPTScc="${optname1} ${optname2} ${optname3} ${optname4}"
OPTScc="${!optname1} ${!optname2} ${!optname3} ${!optname4}"
echo OPTScc=$OPTScc
unset optname1 optname2 optname3 optname4

if test -z "$OPENMP_CFLAGS" ;then
    OPTScc+=" -pthread" #enables pthreads feature and links in the library
    AC_DEFINE(HAS_OPENMP,[0],["Openmp is disabled"])
else
    OPTScc+=" $OPENMP_CFLAGS"
    AC_DEFINE(HAS_OPENMP,[1],["Openmp is enabled"])
fi
OPTSextra=""

#Used for embedding conf folder in maos exe
AC_SUBST(objcopy_output)
AC_SUBST(objcopy_arch)
#AC_CHECK_LIB([fftw3], [fftw_execute],[has_fftw_system="yes"],[has_fftw_system='no'])
#unset ac_cv_lib_fftw3_fftw_execute
#Check operating system.
AM_CONDITIONAL(LINUX,  false)
AM_CONDITIONAL(BSD,    false)
AM_CONDITIONAL(DARWIN, false)
AM_CONDITIONAL(WIN,    false)
mkl_avail="no"
echo host_os=$host_os
case "$host_os" in
    *linux*)
	AM_CONDITIONAL(LINUX, true)
	LDEXE+=" -lrt" #need real time library.
	LDOPTS+=" -Wl,--no-as-needed" #as-needed is default in some system (ubuntu) and cause mkl library checking to fail
	CFLAGS+=" -fPIC"
	CFLAGS_EXTRA="-g -D_XOPEN_SOURCE=500 -D_BSD_SOURCE -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200112L"
	if test "$use_release" = "yes" ;then
	    CFLAGS+=" -D_GLIBC_2_5_SOURCE"
	fi
	if test "${host_cpu}" = "x86_64" ;then
	    libsuffix="64"
	else
	    libsuffix="32"
	fi
	mkl_avail="yes"
	system="linux"
	ldsuffix="so"
	;;
    *darwin*)
	AM_CONDITIONAL(DARWIN, true)
	CFLAGS+=" -fPIC"
	CFLAGS_EXTRA=
	mkl_avail="no"
	libsuffix=fat
	system="apple"
	ldsuffix="dylib"
	#for compatibility with OS X >=10.6. Don't use 10.5 which will render mex files not usable.
	if test "$use_release" = "yes" ;then
	    CFLAGS+=" -mmacosx-version-min=10.7"
	fi
	;;
    *bsd*)
	AM_CONDITIONAL(BSD, true)
	system="bsd"
	ldsuffix="so"
	;;
    *cygwin*|*mingw*)
	AM_CONDITIONAL(WIN, true)
	system="win"
	ldsuffix="dll"
	;;
    none)
	AC_MSG_NOTICE([OS is none])
	;;
    *)
	echo host_os=$host_os
	AC_MSG_ERROR([Unrecognized OS])
esac

AC_DEFINE_UNQUOTED(LDSUFFIX, ["$ldsuffix"], "Suffix of dynamic load library")
PRE_DIR="${HOME}/.aos/"          #prefix for installation
TMP_DIR="${BUILD_DIR}/external/" #temporary directory for compiling
LIB_DIR="${PRE_DIR}/lib${libsuffix}" #external library
LDFLAGS+=" -L${LIB_DIR} -Wl,-rpath,${LIB_DIR}" #-avoid-version
CFLAGS+=" -I${PRE_DIR}/include "
DOWNLOAD_DIR="${TMP_DIR}/downloads"
OPTS+=" $CFLAGS_EXTRA "
if test "$system" = "win" ;then
    LDFLAGS+=" -Wl,-rpath,${PREDIR}/bin"
fi
if test "$host_cpu" = "x86_64" ;then
    CFLAGS+=" -DDLONG "
    mflag="-m64"
else
    CFLAGS+=" -D_FILE_OFFSET_BITS=64" #to fix error about stat() in 32bit.
    mflag="-m32"
fi
if ! echo $CFLAGS | grep -q -- "$mflag" ;then #don't repeat mflag
    CFLAGS=" $mflag $CFLAGS"
fi
unset mflag

die(){
    shift
    echo "Error:$@"
    exit
}


#Discover default download method
if test -x "`which curl 2>/dev/null`" ;then
    wget="curl --retry 1 -L --progress-bar"
    wgetoutput="--output"
elif test -x "`which wget 2>/dev/null`" ;then
    wget="wget -t1 --no-check-certificate"
    wgetoutput="-O"
else
    wget="die wget/curl not found ; "
fi


AM_CONDITIONAL(USE_SCHEDULER, true)
if test "$use_scheduler" = "no" ;then
    AC_DEFINE(MAOS_DISABLE_SCHEDULER, [1], "Disable built in scheduler")
    AM_CONDITIONAL(USE_SCHEDULER, false)
    use_websockets="no"
fi
if test "$host_cpu" = "x86_64" ;then 
    chol_fun=cholmod_l_factorize
    AM_CONDITIONAL(USE_LONG, true)
else
    chol_fun=cholmod_factorize
    AM_CONDITIONAL(USE_LONG, false)
fi

LDEXE+=" -static" #show that we prefer static linking.
if test "x$use_static" = "xyes" ;then
    LDEXE_MAOS="-all-static" #only all static linking for maos
    AM_CONDITIONAL(USE_STATIC, true)
    AC_DEFINE(USE_STATIC, [1], "We are compling maos statically")
else
    AM_CONDITIONAL(USE_STATIC, false)
    AC_DEFINE(USE_STATIC, [0], "We are compling maos dynamically")
fi

if test x$use_oprof = "xyes" ;then
    OPTSprof+=" -g -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls"
fi

#call convention of fortran libraries
#AC_F77_FUNC(dgemm)
#AC_F77_FUNC(dposv)
#download LIBNANE file.tar.bz2 [outdir]
extract(){
    fn="$1"
    OUTDIR="$2"
    if test -z "$OUTDIR" ;then
	AC_MSG_ERROR([Please specify third argument])
	return 1
    else
	echo Extracting "${DOWNLOAD_DIR}/${fn}" to "${OUTDIR}"
	mkdir -p "${OUTDIR}" || die "failed to mkdir"
	if ! tar -xf "${DOWNLOAD_DIR}/${fn}" -C "${OUTDIR}" ;then
	    rm -rf ${DOWNLOAD_DIR}/${fn}
	    AC_MSG_NOTICE([File downloaded is not correct.])
	    return 1
	else
	    return 0
	fi
    fi
}
download(){
    fn="$2"
    OUTDIR="$3"
    mkdir -p ${DOWNLOAD_DIR} || die "failed to mkdir"
    if test -f ${DOWNLOAD_DIR}/${fn} ;then
	extract "${fn}" "$OUTDIR"
    fi
    if ! test -f ${DOWNLOAD_DIR}/${fn} ;then
	echo Downloading $1 from the MAOS website to ${DOWNLOAD_DIR}/${fn}
	for BASEURL in $BASEURLS ; do
	    echo $wget "${BASEURL}/${fn}" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}" 
	    if $wget "${BASEURL}/${fn}" "$wgetoutput" "${DOWNLOAD_DIR}/${fn}" ;then
		if extract "$fn" "$OUTDIR" ;then
		    break
		fi
	    else
		rm -rf ${DOWNLOAD_DIR}/${fn}
		AC_MSG_NOTICE([Failed to download the file from ${BASEURL}/${fn}])
	    fi
	done
	if test ! -f ${DOWNLOAD_DIR}/${fn} ;then
	    AC_MSG_NOTICE(["Failed to download. Please download it manually and put in ${DOWNLOAD_DIR}/${fn}"])
	fi
    fi
    
}
#compile LIBNAME srcfile.tar.gz [compile options]
compile(){
    echo Download and compile $1
    name=$1
    file=$2
    shift
    shift
    rm -rf ${TMP_DIR}/$name
    download $name $file ${TMP_DIR}/$name
    fnlog=${TMP_DIR}/compile.log
    echo "Compiling $name. Check $fnlog for details"
    echo CFLAGS=$CFLAGS > $fnlog
    if (cd ${TMP_DIR}/$name/* && export CFLAGS="$CFLAGS -O3"  && \
	(./configure --libdir=${LIB_DIR} $@ --prefix=${PRE_DIR} >> $fnlog 2>&1|| \
	cmake -DCMAKE_INSTALL_PREFIX:PATH=${PRE_DIR} $@ >> $fnlog 2>&1 ) && \
	make -j4 install >> $fnlog 2>&1) ;then
	echo "Compiled $name successfully."
    else
	echo "Failed to compile $name. Please see $fnlog for details"
    fi
}
function set_lapack(){
    type=$1
    case $type in
	mkl)
	    if test "$type" = "mkl-download" -a "$mkl_avail" = "yes";then
		download MKL "mkl${libsuffix}.tar.bz2" ${PRE_DIR}
	    fi
	    #First use MKL if available. 
	    LAPACK=
	    if test "$host_cpu" = "x86_64" ;then
		LAPACK+=" -lmkl_intel_lp64"
	    else
		LAPACK+=" -lmkl_intel"
	    fi
	    LAPACK+=" -lmkl_core"
	    if test -n "$OPENMP_CFLAGS" ;then
		if test "$CN" = "ICC" ;then
		    LAPACK+=" -lmkl_intel_thread $OPENMP_CFLAGS"
		elif test "$CN" = "GCC" ;then
		    LAPACK+=" -lmkl_gnu_thread $OPENMP_CFLAGS"
		else
		    LAPACK+=" -lmkl_intel_thread $OPENMP_CFLAGS"
		fi
	    else
		LAPACK+=" -lmkl_intel_thread -liomp5"
	    fi
	 
	    if test "$use_static" = "yes" ;then
    		LAPACK="-Wl,--start-group,${LAPACK},--end-group"
		LAPACK=${LAPACK// /,}
	    fi
	    ;;
	pkg-config)
	    LAPACK=`pkg-config --libs blas lapack 2>/dev/null`
	    ;;
	openblas-openmp)
	    LAPACK="-lopenblaso"
	    ;;
	openblas-threads)
	    LAPACK="-lopenblasp"
	    ;;
	openblas)
	    LAPACK="-lopenblas"
	    ;;
	default*)
	    LAPACK="-llapack -lblas "
	    if test "$type" = "default-link" ;then
	       if test -d $/usr/lib${libsuffix} ;then
		   USRLIB=/usr/lib${libsuffix}
	       else
		   USRLIB=/usr/lib
	       fi
	       echo Try to link from $USRLIB
	       fn_blas=`ls ${USRLIB}/libblas.so.* 2>/dev/null`
	       echo fn_blas=${fn_blas}
	       if test -n "${fn_blas}" ;then
		   rm -rf "${LIB_DIR}/libblas.so"
		   ln -sf "${fn_blas%% *}" "${LIB_DIR}/libblas.so"
	       fi
	       fn_lapack=`ls ${USRLIB}/liblapack.so.* 2>/dev/null`
	       if test -n "${fn_lapack}" ;then
		   rm -rf "${LIB_DIR}/liblapack.so"
		   ln -sf "${fn_lapack%% *}" "${LIB_DIR}/liblapack.so"
	       fi
	    fi
	       ;;
		*)
		    LAPACK=
		    ;;
    esac
    echo trying LAPACK=$LAPACK
}
has_lapack=no
base_lapack="pkg-config openblas-openmp openblas-threads openblas default default-link"
case $use_mkl in
    yes) #prefer MKL
	trys="mkl mkl-download $base_lapack"
	;;
    no) #disable MKL
	trys="$base_lapack"
	;;
    default) #default order
	trys="mkl $base_lapack mkl-download"
	;;
esac

for try in ${trys}; do
    if test "$has_lapack" = "no" -o "$has_blas" = "no";then
	unset nset ac_cv_lib_m_dgemm_
	unset nset ac_cv_lib_m_dposv_
	set_lapack ${try}
	if test -n "$LAPACK" ;then
	    #the trailing -lm is necessary for ubuntu-16.04 with mkl
	    AC_CHECK_LIB([m],[dposv_],[has_lapack=yes],[has_lapack=no],[$LAPACK -lm])
	    AC_CHECK_LIB([m],[dgemm_],[has_blas=yes],[has_blas=no],[$LAPACK -lm])
	fi
    fi
done

if test "$has_lapack" = "no" -o "$has_blas" = "no" ;then
    AC_MSG_ERROR([Please install blas/lapack])
fi

#Check for zlib
AC_CHECK_LIB([z],[gzopen],[],[has_z='no'])
if test "$has_z" = "no" ;then
    compile ZLIB zlib.tar.bz2 --static
    unset ac_cv_lib_z_gzopen
    AC_CHECK_LIB([z],[gzopen],[],[AC_MSG_ERROR([Invalid zlib])])
fi

#Check for FFTW interface in MKL
AC_CHECK_LIB([m], [fftw_execute], [builtin_fftw="yes"], [builtin_fftw="no"], [$LAPACK -lm])
if test "$builtin_fftw" = "yes";then
    AC_MSG_NOTICE([FFTW is already provided by $LAPACK])
    fftw_threads=1
else
    #Check for FFTW3 header and library.
    AC_CHECK_LIB([fftw3], [fftw_execute],[has_fftw="yes"],[has_fftw='no'])
    AC_CHECK_HEADERS([fftw3.h],[], [has_fftw='no'])
    if test "$has_fftw" = "no" ;then #fftw does not exist in system.
	FLAGS="--enable-threads --enable-shared --enable-static"
	if test -n "$OPENMP_CFLAGS" ;then
            FLAGS+=" --enable-openmp"
	fi
	compile FFTW fftw.tar.bz2 $FLAGS
	#compile FFTW fftw.tar.bz2 --enable-float 
	unset ac_cv_lib_fftw3_fftw_execute
	AC_CHECK_LIB([fftw3], [fftw_execute],[],[AC_MSG_ERROR([Invalid FFTW library])])
    fi

    if test "$with_matlab" != "no" ;then
	#cannot compile with single precision floating fft when compiling mex routines
	has_fftwf=0
    else
	AC_CHECK_LIB([fftw3f], [fftwf_execute],[has_fftwf=1],[has_fftwf=0])
    fi
    AC_DEFINE_UNQUOTED(HAS_FFTWF, [$has_fftwf], "FFTW single precision exists")
    if test x$has_fftwf = x1 ;then
	LDFFT="-lfftw3 -lfftw3f"
    else
	LDFFT="-lfftw3"
    fi
    fftw_threads=0 #To avoid crush in matlab, we use optional dynamic loading from sys/fft.c instead.
fi
AC_DEFINE_UNQUOTED(BUILTIN_FFTW_THREADS, [$fftw_threads], "Built-in FFTW3 threads support.")

#Check for cholmod
AC_CHECK_LIB([cholmod], [${chol_fun}], [], [has_chol="no"], [$LAPACK -lm])
AC_CHECK_HEADERS([cholmod.h], [], [has_chol="no"])
if test "$has_chol" = "no" ;then
    compile CHOLMOD cholmod.tar.bz2
    unset ac_cv_lib_cholmod_${chol_fun}
    AC_CHECK_LIB([cholmod], [${chol_fun}], [], [AC_MSG_ERROR([Invalid CHOLMOD library])], [$LAPACK -lm])
fi
LDCHOL="-lcholmod"

#Check for libwebsocket
AM_CONDITIONAL(HAS_LWS, false)
if test "$use_websockets" = "yes" ;then
    AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
    AC_CHECK_HEADERS([libwebsockets.h], [], [has_lws="no"])
    if test "$has_lws" = "no" ;then #first try to compile from source. Needs cmake
	compile LIBWEBSOCKETS libwebsockets.tar.bz2 -DLIB_SUFFIX=${libsuffix} -DLWS_WITH_SSL=OFF -DLWS_IPV6=OFF -DLWS_WITHOUT_TESTAPPS=ON -DLWS_WITHOUT_DAEMONIZE=ON -DLWS_WITHOUT_CLIENT=ON -DLWS_WITH_SHARED=OFF -DLWS_WITH_STATIC=ON
	unset ac_cv_lib_websockets_lws_service
	AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
    fi
    if test "$has_lws" = "no" ;then #then download pre-compiled library and header.
	download LIBWEBSOCKETS "libwebsockets${libsuffix}.tar.bz2" ${PRE_DIR}
	unset ac_cv_lib_websockets_lws_service
	AC_CHECK_LIB([websockets], [lws_service], [has_lws="yes"], [has_lws="no"], [])
    fi
    if test "$has_lws" != "no" ;then
	AM_CONDITIONAL(HAS_LWS, true)
	AC_DEFINE(HAS_LWS, [1], "libwebsocket is available")
    else
	AC_MSG_NOTICE(["libwebsocket library or header is not available"])
    fi
fi
AM_CONDITIONAL(DRAW,false) #for drawdaemon
AM_CONDITIONAL(MON, false)  #for parital monitor without notify
AM_CONDITIONAL(MONN,false) #for full monitor with notify
if test "$with_gtk" != "no" ;then
    #Check availability of gtk and libnotify. If available, we will build drawdaemon and the monitor.
    #cairo 1.2 comes with gtk+-2.10
    #make sure there are spaces before and after >=
    if test -n "$GTKDIR" ;then #User specified GTK directory
	if test -d "$GTKDIR/pkgconfig" ;then
	    GTKDIR="$GTKDIR/../"
	fi
	if test -d "$GTKDIR/lib/pkgconfig" ;then
	    export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$GTKDIR/lib/pkgconfig"
	fi
    fi
    PKG_CHECK_MODULES(GTK3, gtk+-3.0,           [gtk_3="yes"], [gtk_3="no"])
    PKG_CHECK_MODULES(GTK2, gtk+-2.0 >= 2.2,    [gtk_2="yes"], [gtk_2="no"])

    #if test x$($PKG_CONFIG --variable=target gdk-2.0) = xquartz ;then
#	macint=gtk-mac-integration
    #fi
    if test "$gtk_2" = "yes" -a "$use_gtk_3" != "yes"; then #prefer gtk+-2
	PKG_CHECK_MODULES(DRAW, gtk+-2.0 >= 2.2 $macint cairo >= 1.4,     [comp_draw="yes"],[comp_draw="no"])
	PKG_CHECK_MODULES(MON,  gtk+-2.0 >= 2.2 $macint ,                 [comp_mon="yes"], [comp_mon="no"])
	PKG_CHECK_MODULES(MONN, gtk+-2.0 >= 2.2 $macint libnotify >= 0.1, [comp_monn="yes"],[comp_monn="no"])
	gtk_ver=2
    elif test "$gtk_3" = "yes" ;then
	PKG_CHECK_MODULES(DRAW, gtk+-3.0 >= 0.0 cairo >= 0.0,     [comp_draw="yes"],[comp_draw="no"])
	PKG_CHECK_MODULES(MON,  gtk+-3.0 >= 0.0,                  [comp_mon="yes"], [comp_mon="no"])
	PKG_CHECK_MODULES(MONN, gtk+-3.0 >= 0.0 libnotify >= 0.0, [comp_monn="yes"],[comp_monn="no"])
	gtk_ver=3
    fi
    if test -n "$GTKDIR" ;then
	TT=-I$GTKDIR/include
	GTK_CFLAGS="$TT $TT/cairo $TT/atk-1.0 $TT/gtk-${gtk_ver}.0 $TT/glib-2.0 $TT/pango-1.0 $TT/gdk-pixbuf-2.0 -I$GTKDIR/lib/glib-2.0/include -I${GTKDIR}/lib/gtk-2.0/include"
	GTK_LIBS="-L$GTKDIR/lib/ -Wl,-rpath,$GTKDIR/lib/"
	if test -d "$GTKDIR/bin" ;then
	    GTK_LIBS="$GTK_LIBS -L$GTKDIR/bin -Wl,-rpath,GTKDIR/bin"
	fi
    else
	GTK_CFLAGS=""
	GTK_LIBS=""
    fi
    if false ;then
	#Temporarily include GTK_CFLAGS in CFLAGS to enable header search.
	CFLAGS_SAVE=$CFLAGS
	CFLAGS="$CFLAGS $GTK_CFLAGS $MON_CFLAGS"
	AC_CHECK_HEADERS([gtk/gtk.h], [], [comp_draw="no"; comp_mon="no"])
	AC_CHECK_HEADERS([libnotify/notify.h], [], [comp_monn="no"])
	CFLAGS=$CFLAGS_SAVE
    fi
    AM_CONDITIONAL(DRAW, test x$comp_draw = xyes)
    AM_CONDITIONAL(MON,  test x$comp_mon  = xyes)

    AM_CONDITIONAL(MONN, test x$comp_monn = xyes)
    if test "$comp_mon" = "yes" ;then
	if ! echo "$MON_LIBS" | grep -q gthread ;then
	    #OLD GLIB has GTHREAD separated
	    PKG_CHECK_MODULES(GTHREAD, gthread-2.0 >= 2.0, [has_gthread="yes"], [has_gthread="no"])
	    if test -n "$GTHREAD_LIBS" ;then
		DRAW_LIBS="$DRAW_LIBS $GTHREAD_LIBS"
		MON_LIBS="$MON_LIBS $GTHREAD_LIBS"
		MONN_LIBS="$MONN_LIBS $GTHREAD_LIBS"
	    fi
	fi
    fi
    #The following will substitute the variables name 
    AC_SUBST(DRAW_CFLAGS)
    AC_SUBST(DRAW_LIBS)
    AC_SUBST(MON_CFLAGS)
    AC_SUBST(MON_LIBS)
    AC_SUBST(MONN_LIBS)
    AC_SUBST(MONN_CFLAGS)
    AC_SUBST(GTK_CFLAGS)
    AC_SUBST(GTK_LIBS)
fi

AC_CHECK_HEADERS([unistd.h])
OPTS+=" $OPTScommon $OPTScc $OPTSprof $OPTSextra"

if test "$with_matlab" != "no" ;then
    if test ${MEXDIR:0:1} = '~' ;then
	MEXDIR=$HOME:${MEXDIR:1}
    fi
    if test -x "${MEXDIR}/mex" ;then
	MEXDIR="${MEXDIR}"
    elif test -x "${MEXDIR}/bin/mex" ;then
	MEXDIR="${MEXDIR}/bin"
    else
	AC_MSG_ERROR([mex not found at $MEXDIR. please reconfigure with --with-matlab=DIR to specify MATLAB location or --without-matlab to disable])
    fi
    MEXEXT="$($MEXDIR/mexext)"
    if test -z "$MEXEXT" ;then
	AC_MSG_ERROR([Unable to determine mex extension. reconfigure with --without-matlab to disable mex.])
    fi
    case $MEXEXT in
   	mexglx)
	    MEXLIB=glnx86
	    MEXOPT="$OPTS -fPIC -m32 -fexceptions -D_FILE_OFFSET_BITS=64"
	    MEXOPT="${MEXOPT} -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
	    mex_cpu=x86
	    ;;
   	mexa64)
	    MEXLIB=glnxa64
	    MEXOPT="$OPTS -fPIC -m64 -fexceptions -fno-omit-frame-pointer"
	    MEXOPT="${MEXOPT} -DDLONG -D_BSD_SOURCE -D_POSIX_C_SOURCE=200112L"
	    mex_cpu=x86_64
	    ;;
	mexmaci)
	    MEXLIB=maci
	    MEXOPT="$OPTS -fPIC -m32 -fexceptions" 
	    mex_cpu=x86
	    ;;
	mexmaci64)
	    MEXLIB=maci64
	    MEXOPT="$OPTS -fPIC -m64 -fexceptions" 
	    MEXOPT="${MEXOPT} -DDLONG"
	    mex_cpu=x86_64
	    ;;
	mexw64)
	    MEXLIB=win64
	    MEXOPT=""
	    mex_cpu=x86_64
	    ;;
	mexw32)
	    MEXLIB=win32
	    MEXOPT=""
	    mex_cpu=x86
	    ;;
	*)
	    AC_MSG_ERROR([Unknown mex extension "$MEXEXT"])
	    ;;
    esac
    if test "$mex_cpu" != "$host_cpu" ;then
	echo mex_cpu=$mex_cpu
	echo host_cpu=$host_cpu
	echo Architecture mismatch, disable mex support.
	mex_lib="no"
    else       
	LDMEX="-L${MEXDIR}/${MEXLIB} -Wl,-rpath,${MEXDIR}/${MEXLIB} "
	AC_CHECK_LIB([mx], [mxGetPr], [mex_lib="yes"], [mex_lib="no"], [$LDMEX])
    fi
    if test x$mex_lib = xyes ;then
    	MEXOPT="${MEXOPT} -I${MEXDIR}/../extern/include -DMATLAB_MEX_FILE $OPENMP_CFLAGS" 
	LIBMEX=" -lmwlapack -lmwblas -lmx -lmex -lmat -lut -lstdc++ -lz -lm -ldl -lfftw3"
	#some MATLAB distribution has libmwfftw3.so, some only has libmwfftw3.so.3
	#AC_CHECK_LIB([mwfftw3], [fftw_execute], [LIBMEX+=" -lmwfftw3"], [LIBMEX+=" -l:libmwfftw3.so.3"], [$LDMEX])
	LDMEX+="-no-fast-install -static -Xcompiler -shared"
	AM_CONDITIONAL(USE_MEX, true)
    else
	with_matlab="no"
	AC_MSG_NOTICE([mex library test failed. Disable mex])
	AM_CONDITIONAL(USE_MEX, false)
    fi
    AC_SUBST(MEXEXT)
    AC_SUBST(LIBMEX)
    AC_SUBST(LDMEX)
    AC_SUBST(MEXOPT)
else
    AM_CONDITIONAL(USE_MEX, false)
fi

if test -z "$CUDADIR" -o "$CUDADIR" = "no" ;then
    with_cuda="no"
fi
if test "$with_cuda" != "no" ;then
    if test -n "$with_cuda" ;then
	user_enable_cuda=yes
    fi
    #if test "$use_static" = "yes"; then
	#AC_MSG_ERROR("all-static conflicts with cuda support. Please disable all-static")
    #fi
    AC_MSG_NOTICE([Enabling CUDA])
    #Check for nvcc executable
    if test ${CUDADIR:0:1} = '~' ;then
	CUDADIR=$HOME:${CUDADIR:1}
    fi
    NVCC="$CUDADIR/bin/nvcc"
    if test ! -x "$NVCC" ;then
	NVCC="$CUDADIR/nvcc"
	CUDADIR="$(dirname $CUDADIR)"
	if test ! -x "$NVCC" ;then
            AC_MSG_NOTICE([nvcc not found])
	    with_cuda="no"
	fi
    fi
    #Check cor cuda library
    if test -f "${CUDADIR}/lib${libsuffix}/libcudart.$ldsuffix" ;then
	CUDA_L="${CUDADIR}/lib${libsuffix}"
    elif test -f "${CUDADIR}/lib/libcudart.$ldsuffix" ;then
	CUDA_L="${CUDADIR}/lib"
    elif test -f "${CUDADIR}/lib/x64/cudart.lib" ;then
    	CUDA_L="${CUDADIR}/lib/x64"
    else
	CUDA_L=""
    fi
    if test -d "$CUDA_L" ;then
	LIBCUDA="-L$CUDA_L -Wl,-rpath,$CUDA_L"
    else
	LIBCUDA=
	AC_MSG_NOTICE([CUDA library not found])
	with_cuda="no"
    fi
    AC_CHECK_LIB([cudart], [cudaSetDevice], [], [AC_MSG_NOTICE([libcudart not found]);with_cuda="no"], [$LIBCUDA])

    #Check for cuda header
    if test -f "${CUDADIR}/include/cuda.h" ;then
	CUDAOPT="-I${CUDADIR}/include"
    else
	AC_CHECK_HEADERS([cuda.h], [] ,[AC_MSG_NOTICE([Header not found]);with_cuda="no"])
	CUDAOPT=
    fi
    if test "$with_cuda" = "no" ;then
	if test "$user_enable_cuda" = "yes" ;then
       	    AC_MSG_ERROR([Specified cuda DIR is not valid])
	else
            AC_MSG_NOTICE(["Cannot compile with cuda. Please use --with-cuda=DIR in configure to specify CUDA location or disable cuda use --without-cuda"])
	fi
    fi
fi

if test "$with_cuda" != "no" ;then
    AC_MSG_NOTICE(["Use cuda in directory $CUDADIR"])
    cudaver="$($NVCC -V |grep release |cut -d ',' -f 2 |cut -d ' ' -f 3)"
    cudaver="${cudaver/.}"
    echo CUDA toolkit is $cudaver
    NVCC+=" $CFLAGS_CPU"
    if test -n "$CCBIN" -a "$with_ccbin" != "no" ;then
	CCBIN="${CCBIN%% *}"
	NVCC+=" -ccbin ${CCBIN%% *}"
    fi
    if test "$cudaver" -lt 40 ;then
	AC_MSG_ERROR([Require at least CUDA toolkit 4.0 ])
    fi
    if test "$cudaver" -ge 40 -a "$cudaver" -lt 90 ;then
	CUDAOPT+=" -gencode arch=compute_20,code=compute_20" #Fermi
    fi
    if test "$cudaver" -ge 42 ;then
	CUDAOPT+=" -gencode arch=compute_30,code=compute_30" #Kepler. Introduced unified memory
    fi
    if test "$cudaver" -ge 65 ;then
        CUDAOPT+=" -gencode arch=compute_50,code=compute_50" #Maxwell
    fi
    if test "$cudaver" -ge 80 ;then
	CUDAOPT+=" -gencode arch=compute_60,code=compute_60" #Pascal. UM transparent migration.
    fi
    if test "$cudaver" -ge 90 ;then
	CUDAOPT+=" -gencode arch=compute_70,code=compute_70" #Volta. UM access counter. ATS over NVlink. Cooperative Groups.
    fi
    if test "$cudaver" -ge 100 ;then
	CUDAOPT+=" -gencode arch=compute_80,code=compute_80" #?
    fi
    if test "$cudaver" -ge 110 ;then
	CUDAOPT+=" -gencode arch=compute_90,code=compute_90" #?
    fi
    if test "$debug" = "yes" ;then
        CUDAOPT+=" -O0"
    else
	CUDAOPT+=" -O3"
    fi
    CUDAOPT+=" -g -lineinfo -DHAVE_CONFIG_H -I.. "
    #CUDAOPT+="-Xptxas=-v -dlcm=cg"
    if test "$host_cpu" = "x86_64" ;then
        CUDAOPT+=" -DDLONG"
    fi
    CCBINOPT="$OPTScommon $OPTSprof $OPTSextra -fPIC $CFLAGS -Wno-unused-parameter"
    CCBINOPT="${CCBINOPT/-fpermissive/}"
    LIBCUDA+=" -lcurand -lcusparse -lcufft -lcublas -lcudart"
    if $CCBIN -v 2>&1 |grep clang ;then
	CCN=CLANG
	LIBCUDA+=" -lc++"
    elif $CCBIN -v 2>&1 |grep icc ;then
	 CCN=ICC
	 LIBCUDA+=" -lstdc++"
    else
	CCN=GCC
	CCBINOPT+=" -Wno-unused-value -D__STRICT_ANSI__ $"
	LIBCUDA+=" -lstdc++"
    fi
    optname=${ccprefix}_${CCN}
    CCBINOPT+=" ${!optname}"
    if test -n "$OPENMP_CFLAGS";then
        CCBINOPT+=" $OPENMP_CFLAGS"
    fi
else
    cudaver="no"
fi
if test "$with_cuda" != "no" ;then
    if test "$cuda_double" = "yes" ;then
	AC_DEFINE(CUDA_DOUBLE, [1], "Use double precision mode in cuda")
	AM_CONDITIONAL(CUDA_DOUBLE, true)
    else
	AC_DEFINE(CUDA_DOUBLE, [0], "Use single precision mode in cuda")
	AM_CONDITIONAL(CUDA_DOUBLE, false)
    fi
    AM_CONDITIONAL(USE_CUDA, true)
    AC_DEFINE(USE_CUDA, [1], "Has Cuda")
    AC_SUBST(CUDAOPT)
    AC_SUBST(CCBINOPT)
    AC_SUBST(LIBCUDA)
    AC_SUBST(NVCC)
else
    AM_CONDITIONAL(USE_CUDA, false)
    AM_CONDITIONAL(CUDA_DOUBLE, false)
    AC_DEFINE(USE_CUDA, [0], "Don't have Cuda")
fi

AC_CONFIG_FILES( [Makefile
	sys/Makefile
	math/Makefile
	lib/Makefile
	tools/Makefile
	maos/Makefile 
	skyc/Makefile
	mex/Makefile
	cuda/Makefile
	])
if test -d "${srcdir}/test" ;then
    AC_CONFIG_FILES([test/Makefile])
fi

LDOPTS+=" -no-fast-install -avoid-version"
LDEXE+=" -lz -lm -ldl"
LDEXE_FULL="$LDFFT $LDCHOL $LAPACK $LDEXE" #when require libaos
echo LIBS=$LIBS
LIBS="" #somehow LIBS contain -fftw3 -lnvidia-smlA
AC_SUBST(ldsuffix)
AC_SUBST(CN)
AC_SUBST(LDAOS)
AC_SUBST(OPTS)
AC_SUBST(LDCHOL)
AC_SUBST(LDFFT)
AC_SUBST(LAPACK)
AC_SUBST(LDEXE)
AC_SUBST(LDEXE_FULL)
AC_SUBST(LDEXE_MAOS)
AC_SUBST(OPTScommon)
AC_SUBST(SRC_DIR)
AC_SUBST(BUILD_DIR)
AC_SUBST(LIBS)
AC_SUBST(LDOPTS)
AC_OUTPUT
echo SRC_DIR=$SRC_DIR
echo BUILD_DIR=$BUILD_DIR
echo host_os=$host_os
echo BIGENDIAN=$BIGENDIAN
echo use_static=$use_static

echo CC=$CC
echo CFLAGS=$CFLAGS
echo LDFLAGS=$LDFLAGS
echo OPTS=$OPTS
echo LDOPTS=$LDOPTS
echo LDEXE_FULL=$LDEXE_FULL
if test "$with_gtk" != "no" ;then
    echo DRAW_LIBS=$DRAW_LIBS
    echo GTK_LIBS=$GTK_LIBS
fi
if test "$with_cuda" != "no" ;then
    echo NVCC=$NVCC
    echo CUDAOPT=$CUDAOPT
    echo LIBCUDA=$LIBCUDA
fi
if test "$with_matlab" != "no" ;then
    echo MEXDIR=$MEXDIR
    echo LIBMEX=$LIBMEX
    echo MEXOPT=$MEXOPT
fi
echo
echo Compiling options:
echo
if test "$with_gtk" != "no" ;then
    echo "Drawdaemon:      $comp_draw"
    echo "Monitor:         $comp_mon (libnotify: $comp_monn)"
    echo "GTK version:     ${gtk_ver}"
fi
if true ;then
    echo "Websocket:       $has_lws"
    echo "CUDA:            $cudaver"
    echo "MATLAB MEX:      $with_matlab"
fi
echo

echo "*************************************************************************"
echo "*** MAOS Documentation can be found at http://lianqiw.github.com/maos ***"
echo "*** MAOS Source can be obtained at http://github.com/lianqiw/maos     ***"
echo "*************************************************************************"
