AC_INIT([RandomFields], 3.5)

CXX=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX
AC_LANG(C++)


CPU_FLAGS=""
AC_CHECK_PROG(LSCPU, [lscpu], yes, no)

if test "x${LSCPU}" == xyes; then
  INFO_CPU=`lscpu | grep Flags | tr "[[:upper:]]" "[[:lower:]]"`
else
  AC_CHECK_PROG(SYSCTL, [sysctl], yes, no)
  if test  "x${SYSCTL}" == xyes; then
    INFO_CPU=`sysctl -a 2> /dev/null | grep machdep.cpu.features | tr "[[:upper:]]" "[[:lower:]]"`
  fi
fi



AC_MSG_CHECKING([which SIMD flags can be recognized])
FX=""
## sysctl -a | grep "cpu.features:" # for OS X
for SET in sse2 ssse3 avx2; do 
  AUTOCPPFLAGS="-m$SET"
  info=`echo "${INFO_CPU}" | grep " $SET "`
  if test "${info}" != ""; then
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]])], [ FX="$FX $AUTOCPPFLAGS" ], [ ])
  fi
done


AUTOCPPFLAGS="$FX"
 
##
AUTOCPPFLAGS="" ## to be deleted if the above does not cause any problems
if test "$AUTOCPPFLAGS" != "" ; then
  AC_MSG_RESULT([$AUTOCPPFLAGS])
else
  AC_MSG_RESULT([none])
fi


RANDOMFIELDS_CXXFLAGS="$CXX_FLAGS $AUTOCPPFLAGS"
RANDOMFIELDS_LIBFLAGS="$LIB_FLAGS"
AC_SUBST(RANDOMFIELDS_CXXFLAGS)
AC_SUBST(RANDOMFIELDS_LIBFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

