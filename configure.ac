# Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

AC_INIT([lobSTR],
	m4_esyscmd([config/git-version-gen .tarball-version]),
        [Melissa Gymrek mgymrek@mit.edu],
        [lobSTR])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([dist-bzip2 -Wall tar-pax foreign subdir-objects])

# Set variables using shell commands
GIT_VERSION=${PACKAGE_VERSION}
AC_SUBST(GIT_VERSION)
MACHTYPE=`uname -m`
AC_SUBST(MACHTYPE)

OS=`uname`
AC_SUBST(OS)
if test "x$OS" = "xDarwin"
then
VARSTATIC=""
AC_SUBST(VARSTATIC)
else
VARSTATIC="-static"
AC_SUBST(VARSTATIC)
fi

CXXFLAGS="$CXXFLAGS"
CFLAGS="$CFLAGS"

AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_FC
AC_PROG_RANLIB
AC_PROG_INSTALL

PKG_CHECK_MODULES([FFTW], [fftw3])
ACX_PTHREAD()
LIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"
AX_CHECK_ZLIB()
AX_BLAS()
if test "$ax_blas_ok" = "no"; then
   AC_MSG_ERROR([
------------------------------------------
Cannot find BLAS libraries.

On Debian/Ubuntu systems, try to install one of the following packages:
   libatlas-dev   liblapack-dev  libblas-dev

If BLAS is installed in a non-standard location,
try adding '--with-blas=<DIR>' (see ./configure --help).
------------------------------------------
])
fi
AX_BOOST_BASE([1.38.0])
PKG_CHECK_MODULES([CPPUNIT],[cppunit])

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
m4_pattern_allow([AM_PROG_AR])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

user_CFLAGS=${CFLAGS}
generic_CFLAGS="-Wall -Wno-strict-aliasing -gdwarf-2 -Wuninitialized"
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

CFLAGS="${generic_CFLAGS} ${ext_CFLAGS} ${user_CFLAGS} ${debug_CFLAGS} ${OPENMP_CFLAGS} -D_GIT_VERSION=\"\\\"${GIT_VERSION}\\\"\" -D_MACHTYPE=\"\\\"${MACHTYPE}\\\"\" -O3"
CXXFLAGS="$CFLAGS"
CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS"
LDFLAGS="$ext_LDFLAGS" 

AC_CONFIG_FILES([
   m4/Makefile
   Makefile
   src/Makefile
   scripts/Makefile
   models/Makefile
])

AC_OUTPUT

dnl
dnl Report Configuration
dnl
AC_MSG_RESULT([])
AC_MSG_RESULT([ Configuration summary for $PACKAGE])
AC_MSG_RESULT([    version $VERSION])
AC_MSG_RESULT([])
AC_MSG_RESULT([ Default installtion directories:])
AC_MSG_RESULT([    programs: ${prefix}/bin/ ])
AC_MSG_RESULT([    support scripts: ${prefix}/share/$PACKAGE/scripts])
AC_MSG_RESULT([    pre-computed model files: ${prefix}/share/$PACKAGE/models])
AC_MSG_RESULT([])
AC_MSG_RESULT([ To change installtion path, use:])
AC_MSG_RESULT([    ./configure --prefix NEW-PATH])
AC_MSG_RESULT([])
