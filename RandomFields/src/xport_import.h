/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 Martin Schlather

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


#ifndef RFxport_H
#define RFxport_H 1

#include "def.h"
#include <zzz_RandomFieldsUtils.h>



#define UTILSCALLS \
  CALL(solve_DELETE);				\
  CALL(solve_NULL);				\
  CALL(solvePosDef);				\
  CALL(solvePosDefSp);				\
  CALL(sqrtPosDefFree);				\
  CALL(sqrtRHS);				\
  CALL(I0mL0);					\
  CALL(invertMatrix);				\
  CALL(getUtilsParam);				\
  CALL(attachRFoptions);			\
  CALL(detachRFoptions);			\
  CALL(ordering);				\
  CALL(orderingInt);				\
  CALL(logWM);					\
  CALL(scalarX);				\
  CALL(detPosDef);				\
  CALL(detPosDefsp);				\
  CALL(XCinvXdet);				\
  CALL(XCinvYdet);				\
  CALL(is_positive_definite);			\
  CALL(chol);					\
  CALL(chol2inv);				\
  CALL(sleepMicro);				\
  CALL(pid);					\
  CALL(utilsparam_DELETE);			\
  CALL(utilsparam_NULL)				\
 
#ifdef CALL
#undef CALL
#endif
#define CALL(what) extern what##_type Ext_##what
UTILSCALLS;

void includeXport();
extern utilsparam* GLOBAL_UTILS;
extern int CORES;
extern int PL;

void PIDKEY_DELETE();
typedef
struct globalparam globalparam;
globalparam *WhichOptionList(bool local);
bool parallel();

//void globalparam_NULL(KEY_type *KT, bool copy_messages);
//void globalparam_NULL(KEY_type *KT);
//void globalparam_DELETE(KEY_type *S);


#endif
