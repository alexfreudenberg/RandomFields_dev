


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#ifndef rfutils_init_H
#define rfutils_init_H 1


 /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
  */

#include "General_utils.h"
#include "Utils.h"
#include "options.h"


#ifdef HAVE_VISIBILITY_ATTRIBUTE
  # define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
  # define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MY_PACKAGE "RandomFieldsUtils"
#define MY_ACRONYM XX
#include "zzz_calls.h"

  /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
!!!!! auch kein mit MALLOC kreiertes Objekt  !!!!
  */

  
  DECLARE1(void, utilsparam_DELETE, utilsparam *, S)
  DECLARE1(void, utilsparam_NULL, utilsparam *, S)
  DECLARE1(void, solve_DELETE, solve_storage**, S)
  DECLARE1(void, solve_NULL, solve_storage*, x)
  DECLARE7(int, solvePosDef, double*, M, int, size, bool, posdef, 
	   double *, rhs, int, rhs_cols, double *, logdet, solve_storage *, PT)
  DECLARE8(int, solvePosDefSp, double *, M, int, size, bool, posdef,
	   double *, rhs, int, rhs_cols, double *,logdet,
	   solve_storage *, Pt, solve_param *,sp)
  //  DECLARE8(int, solvePosDefResult, double*, M, int, size, bool, posdef, 
  //	   double *, rhs, int, rhs_cols, double *, result, double*, logdet, 
  //	   solve_storage*, PT)
  DECLARE4(int, sqrtPosDefFree, double *, M, int, size, solve_storage *, pt,
	   solve_param *, sp)
  DECLARE3(int, sqrtRHS, solve_storage *, pt, double*, RHS, double *, res)
  DECLARE2(int, invertMatrix, double *, M, int, size)
  DECLARE2(double, StruveH, double, x, double, nu)
  DECLARE3(double, StruveL, double, x, double, nu, bool, expScale1d)
  DECLARE1(double, I0mL0, double, x)
  DECLARE3(double, WM, double, x, double, nu, double, factor)
  DECLARE3(double, DWM, double, x, double, nu, double, factor)
  DECLARE3(double, DDWM, double, x, double, nu, double, factor)
  DECLARE3(double, D3WM, double, x, double, nu, double, factor)
  DECLARE3(double, D4WM, double, x, double, nu, double, factor)
  DECLARE4(double, logWM, double, x, double, nu1, double, nu2, double, factor)
  DECLARE1(double, Gauss, double, x)
  DECLARE1(double, DGauss, double, x)
  DECLARE1(double, DDGauss, double, x)
  DECLARE1(double, D3Gauss, double, x)
  DECLARE1(double, D4Gauss, double, x)
  DECLARE1(double, logGauss, double, x)
  
  DECLARE1(void, getUtilsParam, utilsparam **, up)
  DECLARE10(void, attachRFoptions, const char **, prefixlist, int, N, 
	   const char ***, all, int *, allN, setoptions_fctn, set, 
	   finalsetoptions_fctn, final, getoptions_fctn, get,
	    deleteoptions_fctn, del,
	   int, PLoffset,
	   bool, basicopt)
  DECLARE2(void, detachRFoptions, const char **, prefixlist, int, N)

  DECLARE3(void, sorting, double*, data, int, len, usr_bool, NAlast)
  DECLARE3(void, sortingInt, int*, data, int, len, usr_bool, NAlast)
  DECLARE4(void, ordering, double*, data, int, len, int, dim, int *, pos)
  DECLARE4(void, orderingInt, int*, data, int, len, int, dim, int *, pos)
  DECLARE4(double, scalarX, double *, x, double *, y, int, len, int, n)
  //  DECLARE4(int, scalarInt, int *, x, int *, y, int, len, int, n)
  DECLARE2(double, detPosDef, double *, M, int, size) // destroys M!
  DECLARE8(int, XCinvXdet,double*, M, int, size, double *,X, int, X_cols,
	  double *, XCinvX, double *, det, bool, log, solve_storage, *PT)
  DECLARE10(int, XCinvYdet,double*, M, int, size, bool, posdef,
	    double *, X, double *, Y, int, cols,
	    double *, XCinvY, double *, det, bool, log, solve_storage, *PT)
  //  DECLARE5(double, XCinvXlogdet, double *, M, int, size, double *, X,
  //	   int, X_cols, solve_storage *, PT)
  DECLARE2(bool, is_positive_definite, double *, C, int, dim)
  DECLARE2(void, chol2inv, double *, MPT, int, size)
  DECLARE2(int, chol, double *, MPT, int, size)
  DECLARE1(void, pid, int *, i)
  DECLARE1(void, sleepMicro, int *, i)
  // DECLARE7(int, cholGPU, bool, copy, double*, M, int, size, double*, rhs, int, rhs_cols, double *, LogDet, double *, RESULT); // entkommentieren
  

 /* 
!!!!! HIER NIE EIN SEXP OBJEKT ZURUECKGEBEN  !!!!  
  */
 
  /*

    See in R package RandomFields, /src/userinterfaces.cc 
          CALL#(...)
    at the beginning for how to make the functions available
    in a calling package

  */
#ifdef __cplusplus
}
#endif


/*

install.packages("RandomFieldsUtils_0.5.21.tar.gz", configure.args="CXX_FLAGS=-march=native", repo=NULL); library(RandomFieldsUtils, verbose=100); q()

*/

#ifdef DO_PARALLEL
#define HAS_PARALLEL true
#else
#define HAS_PARALLEL false
#endif
#if defined AVX2
#define HAS_AVX2 true
#else
#define HAS_AVX2 false
#endif
#if defined AVX
#define HAS_AVX true
#else
#define HAS_AVX false
#endif
#if defined SSSE3
#define HAS_SSSE3 true
#else
#define HAS_SSSE3 false
#endif
#if defined SSE2
#define HAS_SSE2 true
#else
#define HAS_SSE2 false
#endif
#if defined SSE
#define HAS_SSE true
#else
#define HAS_SSE false
#endif

#define USES_AVX2 (HAS_AVX2 && NEED_AVX2)
#define USES_AVX (HAS_AVX && NEED_AVX)
#define USES_SSSE3 (HAS_SSSE3 && NEED_SSSE3)
#define USES_SSE2 (HAS_SSE2 && NEED_SSE2)
#define USES_SSE (HAS_SSE && NEED_SSE)

#define MISS_AVX2 (!HAS_AVX2 && NEED_AVX2)
#define MISS_AVX (!HAS_AVX && NEED_AVX)
#define MISS_ANY_SIMD (MISS_AVX2 || MISS_AVX || !HAS_SSE2)
#define MISS_SSSE3 (!HAS_SSSE3 && NEED_SSSE3)
#define MISS_SSE2 (!HAS_SSE2 && NEED_SSE2)
#define MISS_SSE (!HAS_SSE && NEED_SSE)


#define AttachMessageN 1000
#define HAS_ONE_RELEVANT (HAS_PARALLEL || USES_AVX2 || USES_AVX || USES_SSSE3 ||  USES_SSE2 || USES_SSE)
#define HAS_ALL_RELEVANT (HAS_PARALLEL && !MISS_AVX2 && !MISS_AVX && !MISS_SSSE3 &&  !MISS_SSE2 && !MISS_SSE)


#define HAS_QUESTIONS							\
  HAS_ONE_RELEVANT ? "sees" : "does not see any of",			\
    HAS_PARALLEL ? "OMP" : "",						\
    USES_AVX2 ? ", AVX2" : "",						\
    USES_AVX ? ", AVX" : "",					\
    USES_SSSE3 ? ", SSSE3" : "",				\
    USES_SSE2 ? ", SSE2" : "",				\
    USES_SSE ? ", SSE" : "",					\
    !HAS_ONE_RELEVANT || HAS_ALL_RELEVANT ? "" : ", but not ",		\
    !HAS_PARALLEL ? "OMP, " : "",					\
    MISS_AVX2 ? "AVX2" : "",				\
    MISS_AVX ? ", AVX" : "",						\
    MISS_SSSE3 ? ", SSSE3" : "",				\
    MISS_SSE2 ? ", SSE2" : "",				\
    MISS_SSE ? ", SSE" : ""


#ifdef WIN32
#define AttachMessageX(PKG, HINT, AND, OMP)				\
  "'"#PKG"' %.20s %.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s.%.320s%.120s%.120s", HAS_QUESTIONS, \
    HINT && MISS_ANY_SIMD ? "\nBy default '"#PKG"' is compiled with flag '-mavx' under your OS.\nIf you are sure that AVX2 is available, consider adding the flag '-march=native'\nto 'PKG_CXXFLAGS' in the file src/Makefile.win and then recompile\n'"#PKG"' "#AND"." : "", \
    HINT && MISS_AVX2 ?							\
    "\nOr: try adding flag '-mavx2' to 'PKG_CXXFLAGS'" : "",\
    HINT && (!HAS_PARALLEL) ? "\nFor OMP alone, try adding the flags -Xpreprocessor -fopenmp -pthread to PKG_LIBS and PKG_CXXFLAGS" : ""
#else
#define AttachMessageX(PKG, HINT, AND, OMP) 				\
  "'"#PKG"' %.20s %.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s%.10s.%.320s%.200s%.200s%.350s%.200s", HAS_QUESTIONS, \
    HINT && MISS_ANY_SIMD ? "\nWithout appropriate SIMD instruction set, the calculations might be slow.\nConsider recompiling '"#PKG"' "#AND" with flags e.g.,\n\n   install.packages(\""#PKG"\", configure.args=\"CXX_FLAGS='-march=native "#OMP"'\")" : "", \
    HINT && MISS_AVX2 ?					\
    "\n\n   install.packages(\""#PKG"\", configure.args=\"CXX_FLAGS='-mavx2 "#OMP"'\")" \
    : "",								\
    HINT && MISS_AVX ?					\
    "\n\n   install.packages(\""#PKG"\", configure.args=\"CXX_FLAGS='-mavx "#OMP"'\")"\
    : "",								\
    HINT && MISS_ANY_SIMD ? "\n\nAlternatively consider installing '"#PKG"'\nfrom https://github.com/schlather/"#PKG", i.e.,\n   install.packages(\"devtools\")\n   library(devtools)\n   devtools::install_github(\"schlather/"#PKG"/pkg\")" : "",\
    HINT && !HAS_PARALLEL ? "\n\nFor OMP alone try\n   install.packages(\""#PKG"\", configure.args=\"CXX_FLAGS='"#OMP"'\")" : ""
#endif

#if defined ownprefixN
#  if defined DO_PARALLEL
#    define AttachMessage(PKG, HINT) AttachMessageX(PKG, HINT, , )
#  else
#    define AttachMessage(PKG, HINT)  \
     AttachMessageX(PKG, HINT, ,\
		    -Xpreprocessor -fopenmp -pthread' LIB_FLAGS='-lomp -pthread)
#  endif
#elif defined DO_PARALLEL
#  define AttachMessage(PKG, HINT)				\
     AttachMessageX(PKG, HINT, and 'RandomFieldsUtils',		\
	            -Xpreprocessor -fopenmp -pthread' LIB_FLAGS='-lomp -pthread)
#else
#  define AttachMessage(PKG, HINT)			\
     AttachMessageX(PKG, HINT, and 'RandomFieldsUtils', )
#endif
  
#define ReturnAttachMessage(PKG,HINT) 	\
  SEXP Ans = PROTECT(allocVector(STRSXP, 1));	\
  char simd_msg[AttachMessageN];			\
  SPRINTF(simd_msg, AttachMessage(PKG,HINT)); \
  SET_STRING_ELT(Ans, 0, mkChar(simd_msg));		\
  UNPROTECT(1);						\
  return Ans;

#endif
