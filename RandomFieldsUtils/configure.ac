
AC_INIT([RandomFieldsUtils], 0.5)
CXX=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX
AC_LANG(C++)



CPU_FLAGS=""
AC_CHECK_PROG(LSCPU, [lscpu], yes, no)

if test "x${LSCPU}" == xyes; then
  INFO_CPU=`lscpu | grep Flags | tr "[[:upper:]]" "[[:lower:]]"`
##   AC_MSG_NOTICE([lspcu yes $INFO_CPU])
 else
  AC_CHECK_PROG(SYSCTL, [sysctl], yes, no)
  if test  "x${SYSCTL}" == xyes; then
##    AC_MSG_NOTICE([lsystl yes])
    INFO_CPU=`sysctl -a 2> /dev/null | grep machdep.cpu.features | tr "[[:upper:]]" "[[:lower:]]"`
  fi
fi



AC_MSG_CHECKING([which SIMD flags can be recognized  cccc])
FX=""
## sysctl -a | grep "cpu.features:" # for OS X
for SET in sse2 ssse3 avx2; do 
  AUTOCPPFLAGS="-m$SET"
##  AC_MSG_NOTICE([flag $AUTOCPPFLAGS])
  info=`echo "${INFO_CPU}" | grep " $SET "`
  if test "${info}" != ""; then
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]])], [ FX="$FX $AUTOCPPFLAGS" ], [ ])
  fi
done

AUTOCPPFLAGS="$FX"
 
##AUTOCPPFLAGS="" ## to be deleted if the above does not cause any problems
if test "$AUTOCPPFLAGS" != "" ; then
  AC_MSG_RESULT([$AUTOCPPFLAGS])
else
  AC_MSG_RESULT([none])
fi

if test "${USE_GPU}" ; then
      CXX_FLAGS="$CXX_FLAGS -DUSEGPU"
      AC_MSG_CHECKING([whether CUDA_HOME is set])
      if test -z "${CUDA_HOME}"; then
          AC_MSG_RESULT(CUDA_HOME not set; using /usr/local/cuda)
          CUDA_HOME="/usr/local/cuda"
      else
          AC_MSG_RESULT(using CUDA_HOME=${CUDA_HOME})
      fi
      AC_CHECK_FILE([${CUDA_HOME}/bin/nvcc], CUDA_INSTALLED=true,)
      if test -z "${CUDA_INSTALLED}"; then
          AC_MSG_RESULT([No CUDA installation found, install CUDA or specify CUDA_HOME.])
      else
          AC_MSG_CHECKING(whether this is the 64 bit linux version of CUDA)
          AC_CHECK_FILE([${CUDA_HOME}/lib64/libcublas.so], , 
          AC_MSG_ERROR([this package only works with 64 bit installations of CUDA]))
          AC_SUBST(CUDA_INSTALLED)

      fi
fi


## install.packages("miraculix_0.9.8.tar.gz", configure.args="CXX_FLAGS=-march=native", repos=NULL)
RANDOMFIELDSUTILS_CXXFLAGS="$CXX_FLAGS $AUTOCPPFLAGS"
RANDOMFIELDSUTILS_LIBFLAGS="$LIB_FLAGS"

#AC_MSG_CHECKING([HuHu])
#AC_MSG_RESULT([$RANDOMFIELDSUTILS_CXXFLAGS])
#AC_MSG_RESULT([$RANDOMFIELDSUTILS_LIBFLAGS])

AC_SUBST(RANDOMFIELDSUTILS_CXXFLAGS)
AC_SUBST(RANDOMFIELDSUTILS_LIBFLAGS)
AC_SUBST(CUDA_HOME)
AC_SUBST(R_LIB)


AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
