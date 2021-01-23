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

//#include "Basic_utils.h"  // must be before anything else
#include "RandomFieldsUtils.h"  // must be before anything else
#ifdef DO_PARALLEL
#include <omp.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include "General_utils.h"
#include "own.h" // vor zzz_RandomFieldsUtils.h
#include "Utils.h"
#include "intrinsics.h"
#include "zzz_RandomFieldsUtils.h" // always last



// local
#ifdef DO_PARALLEL
#else
char ERRMSG[LENERRMSG], MSG[LENERRMSG], MSG2[LENERRMSG];
errorloc_type ERROR_LOC="";
errorstring_type ERRORSTRING;
#endif



int parentpid=-1;
bool parallel() {
  int mypid;
  pid(&mypid);
  return mypid != parentpid;
}

SEXP loadRandomFieldsUtils() {
  pid(&parentpid);
  attachRFoptions(ownprefixlist, ownprefixN,
		  ownall, ownallN,
  		  setparameterUtils,
		  NULL, // final
		  getparameterUtils,
		  delparameterUtils,
		  0, true);
  return R_NilValue;
}


SEXP attachRandomFieldsUtils() {
#define NEED_AVX2 true
#define NEED_AVX true
  //#define NEED_SSE4 true
#define NEED_SSSE3 false
#define NEED_SSE2 true
#define NEED_SSE false
  
  ReturnAttachMessage(RandomFieldsUtils, true);  
}

SEXP detachRandomFieldsUtils(){
  detachRFoptions(ownprefixlist, ownprefixN);
  freeGlobals();
  return R_NilValue;
}



double *ToRealDummy = NULL;
int  ToRealN = 0;
double *ToRealI(SEXP X, bool *create) {
  if (TYPEOF(X) == REALSXP) { 
    *create = false;
    return REAL(X);
  }
  HELPINFO("Better use 'double' as storage mode (for one of the arguments).");
  int len = length(X); 
  double *y;
  if (create || ToRealN < len) {
    y = (double *) MALLOC(sizeof(double) * len);
    if (y == NULL) ERR1("not enough memory for an %d vector of doubles", len);
    if (!create) {
      FREE(ToRealDummy);
      ToRealDummy = y;
      ToRealN = len;
    }
  } else y = ToRealDummy;
  int *x;
  if (TYPEOF(X)==INTSXP) x=INTEGER(X); else x=LOGICAL(X);
  for (int i=0; i<len; i++) y[i] = (double) x[i];
  return y;
}

double *ToReal(SEXP X) {
  bool ToFalse[1] = { false };
 if (TYPEOF(X) == REALSXP) return REAL(X);
  return ToRealI(X, ToFalse);
}


int *ToIntDummy = NULL;
int  ToIntN = 0;
int *ToIntI(SEXP X, bool *create, bool round) {
  if (TYPEOF(X) == INTSXP) {
    *create = false;
    return INTEGER(X);
  }
  if (TYPEOF(X) == LGLSXP) {
    *create = false;
    return LOGICAL(X);
  }
  int len = length(X);
  if (len > 100 || PL > 1)
    HELPINFO("Better use 'integer' as storage mode (for one of the arguments).");
  int *y;
  if (*create || ToIntN < len) {
    y = (int *) MALLOC(sizeof(int) * len);    
    if (y == NULL) ERR1("not enough memory for an %d vector of integers", len);
    if (!*create) {
      FREE(ToIntDummy);
      ToIntDummy = y;
      ToIntN = len;
    }
  } else y = ToIntDummy;
  double *x = (double *) REAL(X);
  if (round) for (int i=0; i<len; i++) y[i] = (int) ROUND(x[i]);
  else for (int i=0; i<len; i++) y[i] = (int) x[i];
  return y;
}

int *ToInt(SEXP X) {
  bool ToFalse[1] = { false };
  return ToIntI(X, ToFalse, false);
}


void freeGlobals() {
  FREE(ToRealDummy);
  FREE(ToIntDummy);
}
