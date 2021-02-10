
AC_INIT([RandomFieldsUtils], 0.6)

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
##CC=`"${R_HOME}/bin/R" CMD config CC`
##CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
##CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
if test -z "$CXX"; then
  AC_MSG_ERROR([No C++ compiler is available])
fi
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
#LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_LANG(C++)
AC_PROG_CXX


AC_MSG_CHECKING([which SIMD flags are available])
##
CPPFLAGS="" ## to be deleted if the above does not cause any problems
if test "$CPPFLAGS" != "" ; then
  AC_MSG_RESULT([$CPPFLAGS])
else
  AC_MSG_RESULT([none])
fi

RANDOMFIELDSUTILS_CXXFLAGS="$CXX_FLAGS $CPPFLAGS"
RANDOMFIELDSUTILS_LIBFLAGS="$LIB_FLAGS"

AC_SUBST(RANDOMFIELDSUTILS_CXXFLAGS)
AC_SUBST(RANDOMFIELDSUTILS_LIBFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

