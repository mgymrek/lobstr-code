# SYNOPSIS
#
#   AX_CHECK_RINSIDE()
#
# DESCRIPTION
#
#   This macro searches for an installed RInside and RCPP library. 

AU_ALIAS([CHECK_RINSIDE], [AX_CHECK_RINSIDE])
AC_DEFUN([AX_CHECK_RINSIDE],

# test if R intsalled 
if test -n "which R"
then
	RCPPFLAGS=$(R CMD config --cppflags)
	RLDFLAGS=$(R CMD config --ldflags)
else
   	AC_MSG_ERROR(R is not installed. Please install R and then continue)
fi

# install Rcpp and RInside
#echo "install.packages('Rcpp',repos='http://cran.r-project.org')" | R --slave --args
#echo "install.packages('RInside',repos='http://cran.r-project.org')" | R --slave --args

RCPPCPPFLAGS=$(echo 'Rcpp:::CxxFlags()' | R --vanilla --slave)
RCPPLDFLAGS=$(echo 'Rcpp:::LdFlags()' | R --vanilla --slave)
RINSIDECPPFLAGS=$(echo 'RInside:::CxxFlags()' | R --vanilla --slave)
RINSIDELDFLAGS=$(echo 'RInside:::LdFlags()' | R --vanilla --slave)
AC_SUBST(RCPPFLAGS)
AC_SUBST(RLDFLAGS)
AC_SUBST(RCPPCPPFLAGS)
AC_SUBST(RCPPLDFLAGS)
AC_SUBST(RINSIDECPPFLAGS)
AC_SUBST(RINSIDELDFLAGS)
)
