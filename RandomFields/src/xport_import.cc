/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by Cholesky or SVD decomposition

 Copyright (C) 2015 -- 2017 Martin Schlather, 

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
#include "RF.h"
#include "xport_import.h"

#define importfrom "RandomFieldsUtils"

#if defined(__clang__)
//# pragma clang diagnostic ignored "-Wcast-function-type"
#endif

#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
#pragma GCC diagnostic ignored "-Wcast-function-type"
#endif


#ifdef CALL
#undef CALL
#endif
#define CALL(what) what##_type Ext_##what = NULL
UTILSCALLS;

#undef CALL
#define CALL(what) Ext_##what = (what##_type) R_GetCCallable(importfrom, #what)
void includeXport() {
  UTILSCALLS;
} // export C

#ifdef __GNUC__
// https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html
#pragma GCC diagnostic warning "-Wcast-function-type"
#endif




int 
PL = 1,
  CORES = INITCORES; //  return;  TO DO: replace by KEYT->global_utils

utilsoption_type *OPTIONS_UTILS;
KEY_type *PIDKEY[PIDMODULUS];
int parentpid=0;
bool parallel() {
  int mypid;
  Ext_pid(&mypid);
  // printf("pid = %d %d\n", mypid, parentpid);
  return mypid != parentpid;
}


void option_NULL(KEY_type *KT, bool keep_messages) {
  //  printf("%d %d %d %d %d; %d %d \n", generalN, gaussN, krigeN, extremeN, fitN, coordsN, prefixN);
  assert(generalN==20 && gaussN == 6 && krigeN == 5 && extremeN == 12
  	 && fitN == 43 && coordsN == 20 && prefixN == 25);
  //assert(MAXDATANAMES == 5 && MAX_COORDNAMES == 4); // increase list in Options.h
  messages_options m;
  if (keep_messages) // unklar wieso dies gebraucht wurde 14.2.21
    MEMCOPY(&m, &(KT->global.messages), sizeof(messages_options));

  MEMCOPY(&(KT->global), &OPTIONS, sizeof(option_type));
  // pointer auf NULL setzten
  if (keep_messages)
    MEMCOPY(&(KT->global.messages), &m, sizeof(messages_options));
    
  Ext_utilsoption_NULL(&(KT->global_utils));
}


void option_DELETE(KEY_type *KT) {
   // pointer loeschen
  Ext_utilsoption_DELETE(&(KT->global_utils));
}



  


void KEY_type_NULL(KEY_type *KT) {
  // JA NICHT UEBER MEMSET, DA PIDs NICHT UEBERSCHRIEBEN WERDEN DUERFEN
 // ******************************************************************
  KT->currentRegister = KT->set = KT->n_data_idx = KT->n_data_names =
    KT->n_coord_idx = KT->n_coord_names = 0;
  KT->naok_range = KT->stored_init = false;
  KT->rawConcerns = unsetConcerns;
  MEMSET(KT->PREF_FAILURE, 0, 90 * Nothing);
  KT->next = NULL;
  KT->error_causing_cov = NULL; // only a pointer, never free it.
  KT->zerox = NULL;
  KT->data_names = KT->coord_names = NULL;
  KT->data_idx =  KT->coord_idx = NULL;
  STRCPY(KT->error_location, "<unkown location>");
  option_NULL(KT, false);
  //   printf("n=%d null=%d\n",  KT->n_coord_names,KT->coord_names!= NULL);
 }

void KEY_type_DELETE(KEY_type **S) {
  KEY_type *KT = *S;
  model **key = KT->KEY;
  option_DELETE(KT);
  FREE(KT->zerox);
    
  for (int nr=0; nr<=MODEL_MAX; nr++)
    if (key[nr]!=NULL) COV_DELETE(key + nr, NULL);

  FREE(KT->data_idx);
  for (int i=0; i<KT->n_data_names; i++) FREE(KT->data_names[i]);
  FREE(KT->data_names);

  FREE(KT->coord_idx);
  for (int i=0; i<KT->n_coord_names; i++) FREE(KT->coord_names[i]);
  FREE(KT->coord_names);
  
  UNCONDFREE(*S);
}


 


KEY_type *KEYT() {  
  int mypid;
  Ext_pid(&mypid);
  KEY_type *p = PIDKEY[mypid % PIDMODULUS];
  if (p == NULL) {
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    assert(neu != NULL);
    assert(neu->zerox == NULL);
    PIDKEY[mypid % PIDMODULUS] = neu;
    neu->visitingpid = mypid;    
    if (PIDKEY[mypid % PIDMODULUS] != neu) { // another process had the
      //                                        same idea
      FREE(neu);
      return KEYT(); // ... and try again
    }
    neu->pid = mypid;
    neu->visitingpid = 0;
    neu->ok = true;
    if (PIDKEY[mypid % PIDMODULUS] != neu) BUG;
    KEY_type_NULL(neu);
    WARN_PARALLEL;
    return neu;
  }
  while (p->pid != mypid && p->next != NULL) p = p->next;
  if (p->pid != mypid) {
    if (!p->ok || p->visitingpid != 0) {
      if (PL >= PL_ERRORS) {
	PRINTF("pid collision %d %d\n",  p->ok, p->visitingpid);
      }
      //    BUG;
     return KEYT();
    }
    p->visitingpid = mypid;
    p->ok = false;
    if (p->visitingpid != mypid || p->ok) return KEYT();
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    neu->currentRegister = UNSET;
    neu->pid = mypid;
    if (!p->ok && p->visitingpid == mypid) {
      p->next = neu;
      p->visitingpid = 0;
      p->ok = true;  
      return neu;
    }
    FREE(neu);
    p->visitingpid = 0;
    p->ok = true;
    KEY_type_NULL(neu);
    return KEYT();
  }
  p->error_causing_cov = NULL;
  return p;
}


SEXP copyoptions() {
  KEY_type *KT = KEYT();
  option_NULL(KT, true);
  return R_NilValue;
}

SEXP setlocalRFutils(SEXP seed, SEXP printlevel) {
  KEY_type *KT = KEYT();
  assert(KT != NULL);
  utilsoption_type *global_utils = &(KT->global_utils);
  assert(global_utils != NULL);
  if (length(seed) > 0)
    global_utils->basic.seed = Integer(seed, (char *) "seed", 0);
  if (length(printlevel) > 0) {
    PL = global_utils->basic.Rprintlevel =
      Integer(printlevel, (char *) "printlevel", 0);
    global_utils->basic.Cprintlevel = global_utils->basic.Rprintlevel +PLoffset;
  }
  return R_NilValue;
}

void finalizeoptions() {
  utilsoption_type *global_utils = OPTIONS_UTILS;
  PL = global_utils->basic.Cprintlevel - PLoffset;
  CORES = global_utils->basic.cores;
}


void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);


#define NEED_GPU false
#define NEED_AVX2 false
#define NEED_AVX true
#define NEED_SSSE3 false
#define NEED_SSE2 true
#define NEED_SSE false
SEXP attachoptions() { // no print commands!!!
#ifdef ReturnAttachMessage
  ReturnAttachMessage(RandomFields, true);
#else
  return R_NilValue;
#endif  
}

void loadoptions() { // no print commands!!!
  for (int i=0; i<PIDMODULUS; i++) PIDKEY[i] = NULL; 
  includeXport();
  Ext_pid(&parentpid);
  Ext_getUtilsParam(&OPTIONS_UTILS);
  utilsoption_type *global_utils = OPTIONS_UTILS;
  global_utils->solve.max_chol = DIRECT_ORIG_MAXVAR;
  global_utils->solve.max_svd = 6555;
  global_utils->solve.pivot_mode = PIVOT_AUTO;
  global_utils->solve.pivot_check = Nan;
  // global_utils->basic.warn_unknown_option = WARN_UNKNOWN_OPTION_NONE1;
  Ext_attachRFoptions((char *) "RandomFields", prefixlist, prefixN, all, allN,
  		      setoptions, finalizeoptions, getoptions, NULL,
  		      PLoffset, true, AVX_NEEDS, GPU_NEEDS, AVX_INFO);
  //  printf("before final\n");
  finalizeoptions();
  //  printf("before init\n"); 
  InitModelList();
}


option_type *WhichOptionList(bool local) {  
  if (local) {
    KEY_type *KT = KEYT();
    if (KT == NULL) BUG;
    return &(KT->global);
  }
  return &OPTIONS;
} 
 
void PIDKEY_DELETE() {
  for (int kn=0; kn<PIDMODULUS; kn++) {
    KEY_type *KT = PIDKEY[kn];
    while (KT != NULL) {
      KEY_type *q = KT;
      KT = KT->next;
      KEY_type_DELETE(&q);
    }
    PIDKEY[kn] = NULL;
  }
}

void detachoptions() {
  PIDKEY_DELETE();
  Ext_detachRFoptions(prefixlist, prefixN);
}


