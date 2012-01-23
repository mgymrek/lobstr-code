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
m4_include([ax_check_rinside.m4])

AC_INIT([lobSTR],
        [1.0.6],
        [Melissa Gymrek mgymrek@mit.edu],
        [lobSTR])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([dist-bzip2])

LF_CONFIGURE_CC
#AC_PROG_CC
#AC_PROG_CPP
#AC_AIX
#AC_ISC_POSIX
AM_PROG_CC_C_O
#AC_MINIX
#AC_HEADER_STDC
LF_CONFIGURE_CXX
#AC_PROG_CXX
#AC_PROG_CXXCPP

LF_HOST_TYPE
LF_SET_WARNINGS
AC_PROG_RANLIB
AC_PROG_CXX

PKG_CHECK_MODULES([FFTW], [fftw3])

ACX_PTHREAD()

LIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"

AX_CHECK_ZLIB()
AX_CHECK_RINSIDE()
AX_BOOST_BASE([1.38.0])



AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_CHECK_FUNCS(fseek64, xpdf_cv_func_fseek64=yes, xpdf_cv_func_fseek64=no)
AC_CHECK_FUNCS(ftell64, xpdf_cv_func_ftell64=yes, xpdf_cv_func_ftell64=no)
if test "$xpdf_cv_func_fseek64" = yes -a "$xpdf_cv_func_ftell64" = yes; then
  AC_DEFINE(HAVE_FSEEK64)
fi



AM_PATH_CPPUNIT(1.12.1)

AC_CONFIG_FILES([
   m4/Makefile
   Makefile
   src/Makefile
])

AC_OUTPUT
