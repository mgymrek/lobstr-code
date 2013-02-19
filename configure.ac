# Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

m4_include([ax_boost_base.m4])
m4_include([ax_check_zlib.m4])
m4_include([ax_blas.m4])

AC_INIT([lobSTR],
        [lobSTR --version or allelotype --version gives current version],
        [Melissa Gymrek mgymrek@mit.edu],
        [lobSTR])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([dist-bzip2])

AM_PATH_CPPUNIT(1.9.6)

# Set variables using shell commands
GIT_VERSION=`git describe --dirty --always`
AC_SUBST(GIT_VERSION)
MACHTYPE=`uname -m`
AC_SUBST(MACHTYPE)

CXXFLAGS="$CXXFLAGS"
CFLAGS="$CFLAGS"

AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_RANLIB
AC_PROG_INSTALL

PKG_CHECK_MODULES([FFTW], [fftw3])
PKG_CHECK_MODULES([PROTOBUF], [protobuf])
ACX_PTHREAD()
LIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"
AX_CHECK_ZLIB()
AX_BLAS()
AX_BOOST_BASE([1.38.0])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_CHECK_FUNCS(fseek64, xpdf_cv_func_fseek64=yes, xpdf_cv_func_fseek64=no)
AC_CHECK_FUNCS(ftell64, xpdf_cv_func_ftell64=yes, xpdf_cv_func_ftell64=no)
if test "$xpdf_cv_func_fseek64" = yes -a "$xpdf_cv_func_ftell64" = yes; then
  AC_DEFINE(HAVE_FSEEK64)
fi

#LT_INIT([disable-shared])

user_CFLAGS=${CFLAGS}
generic_CFLAGS="-Wall -Wno-strict-aliasing -g -gdwarf-2 -Wuninitialized"
ext_CFLAGS=""
debug_CFLAGS=""
case "${host_cpu}-${host_os}" in
  i*86-*linux*)
	    ext_CFLAGS="-march=i686";;
  i*86-darwin*)
	CFLAGS="-m64"
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-arch x86_64"], []);;
  *)
    AC_MSG_CHECKING([if gcc accepts -m64])
    CFLAGS="-m64"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [ext_CFLAGS="-m64"; AC_MSG_RESULT([yes])],
					  [ext_CFLAGS="-D_FILE_OFFSET_BITS=64"; AC_MSG_RESULT([no])]);;
esac

CFLAGS="${generic_CFLAGS} ${ext_CFLAGS} ${user_CFLAGS} ${debug_CFLAGS} ${OPENMP_CFLAGS} -D_GIT_VERSION=\"\\\"${GIT_VERSION}\\\"\" -D_MACHTYPE=\"\\\"${MACHTYPE}\\\"\""
CXXFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS "
LDFLAGS="$ext_LDFLAGS" 

AC_SUBST([VAR_STATIC],["-Wl,-Bstatic"])
AC_SUBST([VAR_DYNAMIC],["-Wl,-Bdynamic"])

AM_INIT_AUTOMAKE([-Wall tar-pax foreign])

AC_CONFIG_FILES([
   m4/Makefile
   Makefile
   src/Makefile
])

AC_OUTPUT
