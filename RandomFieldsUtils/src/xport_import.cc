/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2019 Martin Schlather, 

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
#include "intrinsics.h"
#include "options.h"
#include "xport_import.h"
#include "General_utils.h"
#include "win_linux_aux.h"
#include "zzz_RandomFieldsUtils.h"
#include "RandomFieldsUtils.h"

//#include "zzz_RandomFieldsUtils.h" // always last


#define PLoffset -10

int PL = C_PRINTLEVEL,
  CORES = 1; //  return;  TO DO: replace by KEYT->global_utils

int parentpid=0;
bool parallel() {
  int mypid;
  pid(&mypid);
  return mypid != parentpid;
}

void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);
void deloptions(bool local);


void loadoptions() {
  pid(&parentpid);
  attachRFoptions(prefixlist, prefixN,
		  all, allN,
  		  setoptions,
		  NULL, // final
		  getoptions,
		  deloptions,
		  0, true);
}


SEXP attachoptions() {
#define NEED_AVX2 true
#define NEED_AVX true
  //#define NEED_SSE4 true
#define NEED_SSSE3 false
#define NEED_SSE2 true
#define NEED_SSE false
  ReturnAttachMessage(RandomFieldsUtils, true);  
}



void  detachoptions(){
 freeGlobals();
 detachRFoptions(prefixlist, prefixN);
}


