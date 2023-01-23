#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT([island],[0.2.9])

CC=`"${R_HOME}/bin/R" CMD config CC`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

PKG_CPPFLAGS="-I."

AC_PROG_CC([cc cl gcc])

# GSL configuration based on package 'gsl'
## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
   # Use gsl-config for header and linker arguments
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS=`${GSL_CONFIG} --libs`
else
   AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

CPPFLAGS="${CPPFLAGS} ${GSL_CFLAGS}"
PKG_LIBS="${PKG_LIBS} ${GSL_LIBS}"

AC_MSG_CHECKING([if GSL version >= 1.15])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gsl/gsl_version.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit (major == 1 && minor < 12);
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.15])
else
  AC_MSG_RESULT([yes])
fi

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_CONFIG_HEADERS([src/config.h])
AC_OUTPUT
