
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2017 Martin Schlather, Reinhard Furrer

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

#include <R_ext/Rdynload.h>
#include "def.h"
#include "RandomFields.h"
#include <Basic_utils.h>

#if defined(__clang__)
//# pragma clang diagnostic ignored "-Wcast-function-type"
#endif

#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
// GCC diagnostic ignored "-Wcast-function-type"
#endif

#define none 0
static R_NativePrimitiveArgType
  one_int[] = { INTSXP },
  intdouble[] = {INTSXP, REALSXP }, 
  intintdouble[] = {INTSXP, INTSXP, REALSXP }, 
  charint[] = { STRSXP, INTSXP}, 
  inttwochar[] = {INTSXP, STRSXP, STRSXP }, 
  two_int[] = { INTSXP, INTSXP },
  three_int[] = { INTSXP, INTSXP, INTSXP},
  attr_arg[] = { INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,  INTSXP,
		 INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,  INTSXP };
  //  static R_NativeArgStyle argin[] = {R_ARG_IN},
  //    argout[] = {R_ARG_OUT},
  //   hostarg[] = {R_ARG_OUT, R_ARG_OUT};

#define CDEF(name, n, type) {#name, (DL_FUNC) &name, n, type}
static const R_CMethodDef cMethods[]  = {  
  CDEF(GetAttr, 14, attr_arg),
  CDEF(GetModelName, 3, inttwochar),
  CDEF(GetModelNr, 2, charint),
  CDEF(GetCurrentNrOfModels, 1, one_int),
  CDEF(GetNrParameters, 2, two_int),
  CDEF(PrintModelList, 3, three_int),
  CDEF(PutValuesAtNA, 2, intdouble),
  CDEF(PutValuesAtNAnoInit, 2, intdouble),
  CDEF(expliciteDollarMLE, 2, intdouble),
  CDEF(MultiDimRange, 3, intintdouble),
  CDEF(GetModelRegister, 2, charint),
  //  CDEF(ResetWarnings, 1, one_int),
  CDEF(NoCurrentRegister, 0, none),
  CDEF(GetCurrentRegister, 1, one_int),
  CDEF(PutGlblVar, 2, intdouble),
  CDEF(loadoptions, 0, none),
  CDEF(detachoptions, 0, none),
  {NULL, NULL, 0, none}
};

//SEXP vectordist(SEXP V, SEXP diag);


#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[]  = {
  // in die respectiven C-Dateien muss RandomFieldsUtils.h eingebunden sein
  CALLDEF(attachoptions, 0),
  CALLDEF(copyoptions, 0),
  CALLDEF(DebugCall, 0),
  
  CALLDEF(GetParameterNames, 1),
  CALLDEF(GetSubNames, 1),
  CALLDEF(scatter, 1),
  CALLDEF(GetAllModelNames, 1),
  CALLDEF(GetCoordSystem, 3),
  CALLDEF(GetModelInfo, 5),
  CALLDEF(GetModel, 7),
  CALLDEF(Init, 4),
  CALLDEF(setlocalRFutils, 2),
  CALLDEF(EvaluateModel, 3),
  //  CALLDEF(EvaluateModelXX, 0),
   CALLDEF(GetProcessType, 2),
  CALLDEF(empirical, 10),
  CALLDEF(empvarioXT, 13),
  CALLDEF(fftVario3D, 16),
  CALLDEF(fftVario3DX, 15),
  CALLDEF(boxcounting, 5),
  CALLDEF(detrendedfluc, 5),
  CALLDEF(periodogram, 6),
  CALLDEF(minmax, 5), 
  CALLDEF(MomentsIntern, 2),
  CALLDEF(CovLocNonGrid, 4),
  CALLDEF(LocNonGrid, 2),
  CALLDEF(GetNAPositions, 8),
  CALLDEF(SetAndGetModelFacts, 10),   
  CALLDEF(SetAndGetModelLikelihood, 5), 
  CALLDEF(Take2ndAtNaOf1st, 7),
  CALLDEF(countelements, 3), 
  CALLDEF(countneighbours, 6), 
  CALLDEF(getelements, 5), 
  CALLDEF(getneighbours, 5), 
  CALLDEF(set_boxcox, 2),
  CALLDEF(areDataNamesDifferentFromFormer, 1),
  CALLDEF(areDataIdxDifferentFromFormer, 1),
  CALLDEF(areCoordNamesDifferentFromFormer, 1),
  CALLDEF(areCoordIdxDifferentFromFormer, 1),
  CALLDEF(get_boxcox, 1),
  //  CALLDEF(BoxCox_inverse, 2), 
  CALLDEF(BoxCox_trafo, 4),
  CALLDEF(get_logli_residuals, 1),
  CALLDEF(get_likeliinfo, 1),
  CALLDEF(simple_residuals, 1), 
  CALLDEF(get_linearpart, 2),
  //  CALLDEF(vectordist, 2),
  CALLDEF(maintainers_machine, 0),
//  CALLDEF(),
  {NULL, NULL, 0}
};


void R_init_RandomFields(DllInfo  *dll) {
  R_registerRoutines(dll, cMethods, callMethods, NULL, // .Fortran
		     NULL // extended
		     );
  R_useDynamicSymbols(dll, FALSE); // OK
}


#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
// GCC diagnostic ignored "-Wcast-function-type"
#endif
void R_unload_RandomFields(DllInfo *info) {
  /* Release resources. */
}
#ifdef __GNUC__
// GCC diagnostic warning "-Wcast-function-type"
#endif

