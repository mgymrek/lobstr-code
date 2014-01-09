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
AC_CONFIG_SRCDIR([src/main.cpp])
AM_INIT_AUTOMAKE([dist-bzip2 foreign subdir-objects])

LT_INIT([static])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_FC
AC_PROG_F77

PKG_CHECK_MODULES([FFTW], [fftw3])
AX_PTHREAD()
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

# BamTools needs HAVE_FSEEK64
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_CHECK_FUNCS(fseek64, xpdf_cv_func_fseek64=yes, xpdf_cv_func_fseek64=no)
AC_CHECK_FUNCS(ftell64, xpdf_cv_func_ftell64=yes, xpdf_cv_func_ftell64=no)
if test "$xpdf_cv_func_fseek64" = yes -a "$xpdf_cv_func_ftell64" = yes; then
  AC_DEFINE(HAVE_FSEEK64)
fi

# To compile a static executable (before binary packaging?),
# use:
#    ./configure --enable-all-static
#
# NOTE:
#  1.  Using "./configure LDFLAGS=-static" will not work, because
#      libtool will override the LDFLAG files at times.
#  2.  The flag 'LT_LDFLAGS' is used in ./src/Makefile.am with LDADD.
#
# 'all-static' marco copied from subversion's configure.ac
AC_ARG_ENABLE(all-static,
              AS_HELP_STRING([--enable-all-static],
                             [Build completely static (standalone) binaries.]),
                             [
                              if test "$enableval" = "yes" ; then
                                      LT_LDFLAGS="-all-static $LT_LDFLAGS"
                              elif test "$enableval" != "no" ; then
                                      AC_MSG_ERROR([--enable-all-static doesn't accept argument])
                                      fi
                                      ])

AC_SUBST(LT_LDFLAGS)


# Set variables using shell commands
GIT_VERSION=${PACKAGE_VERSION}
AC_SUBST(GIT_VERSION)
MACHTYPE=`uname -m`
AC_SUBST(MACHTYPE)

# Add pThread related settings,
#  see http://www.gnu.org/software/autoconf-archive/ax_pthread.html
# NOTE: CFLAGS/CXXFLAGS are also updated below.
LIBS="$PTHREAD_LIBS $LIBS"

# Compiler warnings to add
COMPILER_WARNINGS="-Wall -Wextra -Wswitch-default -Wno-strict-aliasing"

# Custom Flags to add to CFLAGS/CXXFLAGS
LOBSTR_CFLAGS="-D_GIT_VERSION=\"\\\"${GIT_VERSION}\\\"\" -D_MACHTYPE=\"\\\"${MACHTYPE}\\\"\""

# Set the final value for CFLAGS/CXXFLAGS
CFLAGS="$CFLAGS $COMPILER_WARNINGS $PTHREAD_CFLAGS $LOBSTR_CFLAGS"
CXXFLAGS="$CXXFLAGS $COMPILER_WARNINGS $PTHREAD_CFLAGS $LOBSTR_CFLAGS $BOOST_CPPFLAGS"

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
