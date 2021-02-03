/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of auxiliary correlation functions 

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is gno error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nonsta     tionary models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc

 Copyright (C) 2005 -- 2017 Martin Schlather

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


 
#include "questions.h"
#include "Processes.h"
#include "families.h"
#include "Coordinate_systems.h"
#include "startGetNset.h"
//#include <R_ext/Lapack.h>
//#include <R_ext/Applic.h>
//#include <R_ext/Utils.h>     
//#include <R_ext/BLAS.h> 




#ifdef SCHLATHERS_MACHINE
#define INFO_TRACE(where, info, cov)
void INFO_TRACEX(const char *where, int *info, model *cov) {
  if (parallel()) BUG;
  static int oldcovnr = UNSET; // ok, da nur schlathers machine
  if (COVNR != oldcovnr && (COVNR < FIRSTDOLLAR || COVNR > LASTDOLLAR)) { 
    oldcovnr = COVNR;
    model *calling = cov->calling;
    if (calling != NULL && (CALLINGNR >= FIRSTDOLLAR && CALLINGNR<=LASTDOLLAR)){
      PRINTF("%s", NAME(calling));
    }
    PRINTF("%s%s (%d i=%d e=%d) [INFO_TRACE]\n", NAME(cov), where,
	   info[INFO_N_X], info[INFO_IDX_X], info[INFO_EXTRA_DATA_X]);
  }
}
#define INFO_TRACE_RETURN(where, info, cov)
void INFO_TRACE_RETURNX(const char *where, int *info, model *cov) {
  if (parallel()) BUG;
  static int oldcovnr = UNSET; // ok, da nur schlathers machine
  if (COVNR != oldcovnr && (COVNR < FIRSTDOLLAR || COVNR > LASTDOLLAR)) { 
    oldcovnr = COVNR;
    model *calling = cov->calling;
    if (calling != NULL && (CALLINGNR >= FIRSTDOLLAR && CALLINGNR<=LASTDOLLAR)){
      PRINTF("%s", NAME(calling));
    }
    PRINTF("RETURN %s%s (%d i=%d e=%d) [INFO_TRACE]\n", NAME(cov), where,
	   info[INFO_N_X], info[INFO_IDX_X], info[INFO_EXTRA_DATA_X]);
  }
}
#else
#define INFO_TRACE(where, info, cov)
#define INFO_TRACE_RETURN(where, info, cov)
#endif



void kdefault(model *cov, int i, double v) {
  utilsparam *global_utils = &(cov->base->global_utils);

  defn *C = DefList + COVNR; // nicht gatternr
  if (PisNULL(i)) {
    switch(C->kappatype[i]) {
    case REALSXP :
      PALLOC(i, 1, 1);
      P(i)[0] = v;
      break;
    case INTSXP :
      PALLOC(i, 1, 1); 
      if (v == NA_INTEGER) PINT(i)[0] = NA_INTEGER;
      else if (!R_FINITE(v)) { BUG }
      else if (v > MAXINT) { BUG}
      else if (v < -MAXINT) { BUG}
      else PINT(i)[0] = (int) v;
      break;
    case STRSXP :
      ERR2("parameter '%.50s' in '%.50s' is undefined.", KNAME(i), NAME(cov));
      break;
    case LISTOF + REALSXP :
      PRINTF("%s:%s (%d) unexpected list\n", NICK(cov), C->kappanames[i], i);
      BUG;
    default : 
      PRINTF("%s:%s (%d) undefined\n", NICK(cov), C->kappanames[i], i);
      BUG;
    }
    cov->nrow[i] = cov->ncol[i] = 1;
  } else if (!global_utils->basic.skipchecks) {    
    if (cov->nrow[i] != 1 || cov->ncol[i] != 1) { 
      LPRINT("%d %s %d nrow=%d, ncol=%d\n", 
	     COVNR, NAME(cov), i, cov->nrow[i], cov->ncol[i]);
      int j; for (j=0; j<cov->ncol[i] * cov->nrow[i]; j++) {
	LPRINT("%10g\n", P(i)[j]);
      }
      ERR3("parameter '%.50s' in '%.50s' is not scalar.%.90s",
	   KNAME(i), NAME(cov), CONTACT);
    }
  }
}


void updatepref(model *cov, model *sub) {
  int i;
  for (i=0; i< Forbidden; i++) {
    if (i==Specific) continue;
    if (sub->pref[i] < cov->pref[i]) {
      cov->pref[i] = sub->pref[i];
    }
  }
}


void setdefault(model *cov, int vdim0, int vdim1) {
  // Achilles-Ferse: setdefault wird aufgerufen bevor die Varianten
  // der Modelle und Operatoren feststehen -- die Parameter unten muessen, falls
  // notwendig von Hand in den checks gesetzt werden.

  defn *C = DefList + COVNR;// nicht gatternr 
  system_type *def = DEF;
  int
    n = SYSTEMS(def),
    last = LASTSYSTEM(def);

  //  printf("'%s' ----------------\n", NAME(cov)); 
  
  cov->full_derivs = C->F_derivs;
  cov->rese_derivs = C->RS_derivs;
  cov->loggiven = (ext_bool) (C->log != ErrLogCov); 
 
  for (int s=0; s<n; s++) {
    set_last(OWN, s, last);
    set_maxdim(OWN, s, MAXDIM(def, s)); 
  }

  cov->logspeed = RF_NA;
  
  if ((C->vdim == PREVMODEL_DEP  ||  C->vdim == SUBMODEL_DEP)) {    
    VDIM0 = vdim0;
    VDIM1 = vdim1;
  } else {
    assert(C->vdim != MISMATCH); // ???
    VDIM0 = VDIM1 = C->vdim;     
  }

  if (isPosDef(cov)) { // kein isnow!!, da dies gebraucht wird, wenn
    // posdef fkten als shape fcts verwendet werden.
    for (int i=0; i<MAXMPPVDIM; i++) cov->mpp.maxheights[i] = 1.0; // maxv
  }

  if (isIsotropic(def) && isIsotropic(OWN) && 
      isPosDef(OWNTYPE(0)) && isXonly(C->systems[0]))
    cov->logspeed = 0.0;
  
  cov->finiterange = C->finiterange;  
  cov->monotone=C->Monotone;
  cov->ptwise_definite=C->ptwise_definite;

  //   cov->diag = 
  //   cov->semiseparatelast = 
  //   cov->separatelast = isIsotropic(OWNISO(0)); // war false vor 27.9.14

  
  MEMCOPY(cov->pref, C->pref, sizeof(pref_shorttype));
  
  //for (Methods i=Nothing+1; i<Forbidden; i++) cov->pref[i] = PREF_NONE;

  cov->method = Forbidden; 

  cov->taylorN = C->TaylorN;
  cov->tailN = C->TailN;

  for (int i=0; i<cov->taylorN; i++) {
    cov->taylor[i][TaylorConst] = C->Taylor[i][TaylorConst];
    cov->taylor[i][TaylorPow] = C->Taylor[i][TaylorPow];
  }
  for (int i=0; i<cov->tailN; i++) {
    cov->tail[i][TaylorConst] = C->Tail[i][TaylorConst];
    cov->tail[i][TaylorPow] = C->Tail[i][TaylorPow];
    cov->tail[i][TaylorExpConst] = C->Tail[i][TaylorExpConst];
    cov->tail[i][TaylorExpPow] = C->Tail[i][TaylorExpPow];
  }
}


//int merge_integer(int val, int sub){
//  if (val == SUBMODEL_DEP) return sub; // for Dimension
//  else {    
//    //    if (*val == PARAM_DEP) {model *cov; crash(cov);}
//    assert(val != PARAM_DEP);
//    if (sub < val) return sub;
//  }
//  return val;
//}


ext_bool merge_extbool(ext_bool val, ext_bool sub){
  if (val == SUBMODEL_DEP) return sub;
  else {
    assert((int) Paramdep < 0 && (int) Submodeldep < 0);
    if ((int) sub < (int) val) return sub;
  }
  return val;
}


monotone_type merge_monotone(monotone_type val, monotone_type sub){
  if (val == MON_SUB_DEP) return sub;
  else {
    if ((int) sub < (int) val) return sub;
  }
  return val;
}


// no setbackard
void setbackward(model *cov, model *sub) {
   defn *C = DefList + COVNR;// nicht gatternr
 // see also check_co when changing

  assert(cov != NULL);
  //if (sub == NULL) crash();
  assert(sub != NULL);
 
  int last = OWNLASTSYSTEM;
  if (last == SUBLASTSYSTEM) {
    for (int s=0; s<=last; s++) { 
      //    set_maxdim(OWN, s, merge_integer(MAXDIM(OWN, s), MAXDIM(SUB, s)));
    }
  }
  cov->monotone = merge_monotone(cov->monotone, sub->monotone);  
  cov->finiterange = merge_extbool(cov->finiterange, sub->finiterange);
  
  if (sub->full_derivs < cov->full_derivs || cov->full_derivs == UNSET)
      cov->full_derivs = sub->full_derivs;

  assert(cov->full_derivs >= 0 || 
	 (cov->full_derivs == MISMATCH && isRandom(cov) && isRandom(sub)));
  if (sub->rese_derivs < cov->rese_derivs || cov->rese_derivs == UNSET)
    cov->rese_derivs = sub->rese_derivs;
  if (cov->loggiven != falsch) cov->loggiven = sub->loggiven;
  
  updatepref(cov, sub);

  if (sub==cov->sub[0] || sub==cov->key) {
    if (C->vdim == SUBMODEL_DEP) {
      VDIM0 = sub->vdim[0];
      VDIM1 = sub->vdim[1];
    }
    if (C->ptwise_definite == pt_submodeldep) {
      cov->ptwise_definite = sub->ptwise_definite;
      // assert((cov->ptwise_definite != pt_paramdep &&  // to do !
      //	      cov->ptwise_definite != pt_submodeldep &&
      //	      cov->ptwise_definite != pt_undefined));
    } 
  } else {
    if (cov->ptwise_definite != sub->ptwise_definite) {
       // assert((cov->ptwise_definite != pt_paramdep && // to do !
       //     cov->ptwise_definite != pt_submodeldep &&
       ///     cov->ptwise_definite != pt_undefined &&
       //     sub->ptwise_definite != pt_paramdep &&
       //     sub->ptwise_definite != pt_submodeldep &&
       //     sub->ptwise_definite != pt_undefined));
      if (cov->ptwise_definite ==pt_mismatch ||
	  sub->ptwise_definite==pt_mismatch)
	cov->ptwise_definite = pt_mismatch;
      else if (cov->ptwise_definite==pt_unknown ||
	       sub->ptwise_definite==pt_unknown)
	cov->ptwise_definite = pt_unknown;
      else 	  
	cov->ptwise_definite = 
	  cov->ptwise_definite == pt_zero ? sub->ptwise_definite
	  : sub->ptwise_definite == pt_zero ? cov->ptwise_definite      
	  : pt_indef;
    }
  }

  cov->hess = (DefList[COVNR].hess != NULL && sub->hess);
  cov->randomkappa |= sub->randomkappa;
}


int alloc_mpp_M(model *cov, int moments) {
  int maxmoments = DefList[COVNR].maxmoments;
  // insbesondere fuer models die selbst vom Random-Type sind
  // printf("alloc %s %d max=%d==%d\n", NAME(cov), moments, maxmoments, PARAM_DEP);
  assert(moments >= 0);
  //  printf("maxmoments %s = %d %d %d\n", NAME(cov), maxmoments, MISMATCH, PARAM_DEP);
  //PMI(cov);
  //if (maxmoments == PARAM_DEP) crash();
  // PMI0(cov);
  assert(maxmoments != MISMATCH && maxmoments != PARAM_DEP);


  if (moments > maxmoments && maxmoments != SUBMODEL_DEP) {
    SERR2("required moments (%d) exceeds the coded moments (%d)",
	  moments, maxmoments);
  }

  
  if (moments <= cov->mpp.moments) RETURN_NOERROR;
  if (cov->mpp.mM != NULL) free_mpp_M(cov);
  cov->mpp.moments = moments;

  int 
    vdim = VDIM0,
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    bytes = sizeof(double) * nmvdim;

  if (vdim <= 0) BUG;
  // if (vdim > MAXMPPVDIM) SERR1("multivariate dimension (%d) too large", vdim);
 
  cov->mpp.mM = (double*) MALLOC(bytes);
  cov->mpp.mMplus = (double*) MALLOC(bytes);

  //  assert(nm < 100);
  int nmP1 = cov->mpp.moments + 1;
  for (int i=0; i<nmvdim; i++) cov->mpp.mMplus[i] = cov->mpp.mM[i] = RF_NA;
  double maxheight = isPosDef(cov) ? 1.0 : RF_NAN;
  for (int v=0; v<vdim; v++) {
    int idx = v * nmP1;
    cov->mpp.mMplus[idx + 0] = cov->mpp.mM[idx + 0] = RF_INF;
    if (v < MAXMPPVDIM) cov->mpp.maxheights[v] = maxheight; // maxv
  }

  // cov->mpp.mMplus[0] = cov->mpp.mM[0] = 1.0;
  RETURN_NOERROR;
}

void free_mpp_M(model *cov) {
  FREE(cov->mpp.mM);
  FREE(cov->mpp.mMplus);
  cov->mpp.mM = cov->mpp.mMplus = NULL;
}

int UpdateMPPprev(model * cov, int moments) {
  model *calling = cov->calling;
  int i, nm, err,
    nmvdim,
    vdim = VDIM0;

  nm = cov->mpp.moments < calling->mpp.moments ? cov->mpp.moments
	: calling->mpp.moments;
  nmvdim = (nm + 1) * vdim;
  if (moments >= 0 && calling != NULL) {
    if (calling->mpp.moments == SUBMODEL_DEP &&
	(err = alloc_mpp_M(calling, moments)) != NOERROR) RETURN_ERR(err);
    for (i=0; i<nmvdim; i++) {
      calling->mpp.mMplus[i] = cov->mpp.mMplus[i];
      calling->mpp.mM[i]     = cov->mpp.mM[i];
    }
  }
  
  // nachfolgende Zeilen so lassen, da sonst unerwuenscht
  // maxheight etc. nach oben gegeben werden.
  // calling->mpp.maxheight = cov->mpp.maxheight; 
  // calling->mpp.unnormedmass = cov->mpp.unnormedmass;
  
  RETURN_NOERROR;
}


int INIT_intern(model *cov, int moments, gen_storage *s) { // kein err

  //  printf("initialising %s\n", NAME(cov));
  
  if (!cov->checked) BUG;
  if (cov->initialised) RETURN_NOERROR;
  assert(cov != NULL);
  ASSERT_GATTER(cov);
  

  defn *C = DefList + COVNR;
  int err = NOERROR;
  char *error_location = cov->base->error_location;  

  SPRINTF(error_location, "initializing %.50s", NICK(cov));

  // printf("Iintern %s = %d mom=%d, vdim=%d\n", NAME(cov), cov->mpp.moments, moments, VDIM0);
  
  if (cov->mpp.moments != SUBMODEL_DEP && cov->mpp.moments != PARAM_DEP) {
    if ((err = alloc_mpp_M(cov, moments)) != NOERROR) RETURN_ERR(err);
  } else {
    BUG; //assert(false); // passt das hier??
    if (cov->mpp.moments == PARAM_DEP) cov->mpp.moments = moments;
  }
 
  if (C->maxmoments >= 0 && moments > C->maxmoments) {
      SERR3("moments known up to order %d for '%.50s', but order %d required",
	    C->maxmoments, NICK(cov), moments);
  }
  
  SPRINTF(error_location, "%.50s", cov->calling == NULL ? "initiating the model"
	  : NICK(cov->calling));
  ASSERT_GATTER(cov);
  if ((err = DefList[GATTERNR].Init(cov, s)) != NOERROR) {
    RETURN_ERR(err);   
  }
   
  if ((err = UpdateMPPprev(cov, moments)) != NOERROR) {
    RETURN_ERR(err);
  }

  cov->initialised = true;
  // printf("returning from initialising %s without e rror.\n", NAME(cov)); // PMI0(cov);

 RETURN_NOERROR;
}

void set_initialised_false(model *cov){
  int i;
  if (!cov->randomkappa) return;
  cov->initialised = false;

  for (i=0; i<MAXPARAM; i++) {
    if (cov->kappasub[i] != NULL) {
      set_initialised_false(cov->kappasub[i]);
    }
  }

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL)
      set_initialised_false(cov->sub[i]);
  }
}


int REINIT_intern(model *cov, int moments, gen_storage *s) { // kein err
  int err;
  set_initialised_false(cov);
  err = INIT_intern(cov, moments, s);
  RETURN_ERR(err);
}


int INIT_RANDOM_intern(model *cov, int moments, gen_storage *s, // kein err
		       double *p) {
  if (!cov->checked) BUG;
  if (!cov->initialised) {
    int err = NOERROR;
    char *error_location = cov->base->error_location;  

    SPRINTF(error_location, "initializing %.50s", NICK(cov));
     
    assert(cov != NULL);
    if (moments < 0) SERR("moments expected to be positive");
    if (DefList[COVNR].maxmoments >= 0 &&
	moments > DefList[COVNR].maxmoments) SERR("Moments do not match");
    
    if (cov->mpp.moments != SUBMODEL_DEP && cov->mpp.moments != PARAM_DEP) {
      if ((err = alloc_mpp_M(cov, moments)) != NOERROR) RETURN_ERR(err);
    } else {
      BUG; // passt das hier??
      if (cov->mpp.moments == PARAM_DEP) cov->mpp.moments = moments;
    }
    
    SPRINTF(error_location, "%.50s", cov->calling == NULL
	    ? "initiating the model" : NICK(cov->calling));
    ASSERT_GATTER(cov);
    if ((err = DefList[GATTERNR].Init(cov, s)) != NOERROR) RETURN_ERR(err);   
    if (ISNAN(cov->mpp.mM[moments])) {
      SERR1("%.50s is not a random function", NICK(cov));
    }
    
    if ((err = UpdateMPPprev(cov, moments)) != NOERROR) RETURN_ERR(err);
    
    cov->initialised = true;
  }

  //  switch (DefList[COVNR].kappatype[param_nr]) {
  //  case REALSXP :   
  if (s->dosimulate) DORANDOM(cov, p);
    

    //    break;
    //  case INTSXP :
    //    int j, len;
    //    double *dummy;
    //    dummy = (double*) MALLOC(sizeof(double) * len);
    //    DORANDOM(cov, dummy);
    //    for (j=0; j<len; j++) p[j] = (int) dummy[j];
    //    FREE(dummy);
    //    break;
    //  default : SERR("random parameter only allowed for numerical values");
    //  }
  
  RETURN_NOERROR;
}


void stat2_Intern(double *x, model *cov, double **Z) {
  //    PMI(cov->calling);
  double *z = *Z;
  TALLOC_DOUBLE(z1);
  int trafonr = TRAFONR;
  bool trafo = cov->calling != NULL && trafonr != UNSET;
  DEFAULT_INFO(info);

   if (trafo) {
    TALLOC_GATTER_GLOBAL(z1, GATTERTOTALXDIM);
    DefList[trafonr].cov(x, info, cov, z1); // Trafo-FCTN() call
    x = z1;
  }
 
  int n = GATTERLASTSYSTEM;
  //  PMI0(cov);  printf("n=%d\n", n);
  assert(n == 0); // **Z is a workaround for more complicated systems (i_nrow)
  for (int s=0; s<=n; s++) {
    int gnr = NRi(GATTER[s]), // OK
      dim = GATTERXDIM(s);
     //    printf("stat2 %s\n", NAME(cov));
    assert(dim > 0 && gnr >= FIRSTGATTER && gnr <= LASTGATTER);
    
    switch(gnr) {
    case ISO2ISO :
      *x=FABS(*x);
      if (trafo) *z = *x; else *Z = x;
      break; // crucuial for i_row, r_col
    case S2ISO : case SP2ISO : {
      double b = 0.0;
      for (int d=0; d<dim; d++) b += x[d] * x[d];
      *z = SQRT(b);
      break;
    }
    case S2S : case SId :
      //      printf("dim = %d %s %d\n", dim, NAME(cov), trafo);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
     break;
    case SP2SP :
      x[0]=FABS(x[0]); x[1]=FABS(x[1]);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
      break;
    case S2SP : {
      double b = 0.0;
      int dimM1 = dim - 1;
      for (int d=0; d<dimM1; d++) b += x[d] * x[d];
      z[0] = SQRT(b);
      z[1] = FABS(x[dimM1]);
      break;
    }      
    case E2EIso : EarthIso2EarthIso(x);
      if (trafo) *z = *x; else *Z = x;
      break;
    case E2E : Earth2Earth(x);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
      break;
    case E2SphIso : EarthIso2SphereIso(x, cov, z); break;
    case E2Sph : Earth2Sphere(x, cov, z); break;
    case Sph2SphIso : SphereIso2SphereIso(x);
      if (trafo) *z = *x; else *Z = x;
      break;    
    case Sph2Sph : Sphere2Sphere(x);
      if (trafo) MEMCOPY(z, x, sizeof(double) * dim); else *Z = x;
      break;  
    default :
      //      printf("gnr=%d\n", gnr);
      //      PMI(cov);
      BUG;
    }
    z += OWNXDIM(s);
    x += dim;
  }
  FREE_TALLOC(z1);
}

//


void stat2(double *x, int *info, model *cov, double *v) {
  //  if (!equalsXonly(PREVDOM(0))) { PMI0(cov); crash(); }
  assert(equalsXonly(PREVDOM(0)));
 // printf("stat2\n");
  INFO_TRACE("2", info, cov);
  TALLOC_GATTER(z, OWNTOTALXDIM);
  double **z1 = &z;
  stat2_Intern(x, cov, z1);
  //  if (OWNTOTALXDIM == 2)
  //printf("z after =%ld %10g %10g %10g; x=%ld %d xdim=%d\n", (*z1),(*z1)[0],(*z1)[1],
  //	   (*z1)[2], x, VDIM0, OWNTOTALXDIM);
  //  printf("stat2: %s\n", NAME(cov));

  
  DefList[COVNR].cov(*z1, info, cov, v);// nicht gatternr
  END_TALLOC_z; // very crucial that not z is freed!
  INFO_TRACE_RETURN("2", info, cov);
}

void logstat2(double *x, int *info, model *cov, double *v, double *Sign) {
  assert(equalsXonly(PREVDOM(0)));
  INFO_TRACE("2L", info, cov);
  TALLOC_GATTER(z, OWNTOTALXDIM);
  double **z1 = &z;
  stat2_Intern(x, cov, z1);
  DefList[COVNR].log(*z1, info, cov, v, Sign);// nicht gatternr
  END_TALLOC_z;
  INFO_TRACE_RETURN("2L", info, cov);
}


void nonstat2stat(double *x, double *y, model *cov, double *z) {
  assert(PREVDOM(0) == KERNEL);
  // printf("nonstat2stat\n");
  int n = GATTERLASTSYSTEM;
  for (int s=0; s<=n; s++) {
    int gnr = NRi(GATTER[s]), // OK
      dim = GATTERXDIM(s);

    //    printf("nonstat2 %s\n", NAME(cov));
    assert(dim > 0 && gnr >= FIRSTGATTER && gnr <= LASTGATTER);
    switch(gnr) {
    case S2ISO : {
      double b = 0.0;
      for (int d=0; d<dim; d++) {
	//	printf("d=%d  %10g %10g %d\n", d, x[d], y[d], dim);
	//	if (dim == 1) APMI(cov->calling);
	
	double a = x[d] - y[d];
	b += a * a;
      }
      *z = SQRT(b);
      break;
    }
    case S2S :case SId:
      // printf("d=%d ", dim);
      for (int d=0; d<dim; d++) {
	//	printf("%d\n", d);
	//	printf("%d x=%f; ", d, x[d]);
	//	printf("%d y=%f; ", d, y[d]);
	z[d] = x[d] - y[d];
      }
      break;
    case S2SP :{
      double b = 0.0;
      int dimM1 = dim - 1;
      for (int d=0; d<dimM1; d++) {
	double a = x[d] - y[d];
	b += a * a;
      }
      z[0] = SQRT(b);
      z[1] = FABS(x[dimM1] - y[dimM1]);
      break;
    }    
    case E2EIso : nonstatEarth2EarthIso(x, y, cov, z); break;
      //    case E2E : nonstatEarth2Earth(x, y, cov, z); break;
    case E2SphIso : nonstatEarth2SphereIso(x, y, cov, z); break;
      //    case E2Sph : nonstatEarth2Sphere(x, y, cov, z); break;
    case Sph2SphIso : nonstatSphere2SphereIso(x, y, cov, z); break;     
      //    case Sph2Sph : nonstatSphere2Sphere(x, y, cov, z); break;  
    default : //e.g. SId, E2E, ...
      // S2FORBIDDEN :
      //      printf("gnr=%d\n", gnr);
      //      PMI(cov);
      BUG;
    }
    z += OWNXDIM(s);
    x += dim;
    y += dim;
  }
}

void nonstat2(double *x, double *y, int *info, model *cov, double *v) {
  //  if (PREVDOM(0) != KERNEL) {PMI0(cov); crash(); }
  assert(PREVDOM(0) == KERNEL);
  //  printf("nonsat2 !!\n");
  INFO_TRACE("N2", info, cov);
  //  *v = 1.0; return;
  TALLOC_DOUBLE(z1);
  TALLOC_DOUBLE(z2);
  int nr = TRAFONR; 
  if (cov->calling != NULL && nr != UNSET) {
    //    printf("earth trafo:\n");
    int
      prevdim = PREVXDIM(0),
      xdim = GATTERTOTALXDIM;
    TALLOC_GATTER_GLOBAL(z1, xdim);
    TALLOC_GATTER_GLOBAL(z2, xdim);
    DefList[nr].cov(x, info, cov, z1);  // Trafo-FCTN() call 
    DefList[nr].cov(y, info, cov, z2); // Trafo-FCTN() call 
    x = z1;
    y = z2;
  }
  if (equalsKernel(OWNDOM(0))) {
    assert(DefList[COVNR].nonstat_cov != nonstat2);
    DefList[COVNR].nonstat_cov(x, y, info, cov, v);// nicht gatternr
  } else {
    // printf("> %s ", NAME(cov));
    // kernel2xonly:
    TALLOC_GATTER(z, OWNTOTALXDIM);
    nonstat2stat(x, y, cov, z);
    DefList[COVNR].cov(z, info, cov, v);// nicht gatternr 
    END_TALLOC_z;
  }
  FREE_TALLOC(z1);
  FREE_TALLOC(z2);
  INFO_TRACE_RETURN("N2", info, cov);
}

void nonstat_log2(double *x, double *y, int *info, model *cov,
		       double *v, double *Sign) {
  assert(PREVDOM(0) == KERNEL);
  INFO_TRACE("NL2", info, cov);
  TALLOC_DOUBLE(z1);
  TALLOC_DOUBLE(z2);
  int nr = TRAFONR;
  if (cov->calling != NULL && nr != UNSET) {
    int xdim = GATTERTOTALXDIM;
    TALLOC_GATTER_GLOBAL(z1, xdim);
    TALLOC_GATTER_GLOBAL(z2, xdim);
    DefList[nr].cov(x, info, cov, z1); // Trafo-FCTN() call 
    DefList[nr].cov(y, info, cov, z2); // Trafo-FCTN() call 
    x = z1;
    y = z2;
  }
  if (equalsKernel(OWNDOM(0))) {
    DefList[COVNR].nonstatlog(x, y, info, cov, v, Sign);// nicht gatternr
    return;
  } else {
    // kernel2xonly:
    TALLOC_GATTER(z, OWNTOTALXDIM);
    nonstat2stat(x, y, cov, z);
    DefList[COVNR].log(z, info, cov, v, Sign);// nicht gatternr
    END_TALLOC_z;
  }
  FREE_TALLOC(z1);
  FREE_TALLOC(z2);
  INFO_TRACE_RETURN("NL2", info, cov);
}


void D_2(double *x, int *info, model *cov, double *v){  
  INFO_TRACE("D", info, cov);
  assert(everyCoord(isSpaceIsotropic, cov));
  defn *C = DefList + COVNR;// nicht gatternr
  int dim = GATTERXDIM(0); // frueher prevxdim
  if (dim == 1) {
    assert(equalsIsotropic(OWNISO(0)));
    double y = FABS(*x);

    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht

    C->D(&y, info, cov, v);// nicht gatternr
  } else {
    assert(dim == 2);
    if (OWNTOTALXDIM == 1) {
      assert(equalsIsotropic(OWNISO(0)));
      double r=SQRT(x[0] * x[0] + x[1] * x[1]);    
      C->D(&r, info, cov, v);
      if (r!=0.0) *v *= x[0] / r;
    } else {
      assert(OWNTOTALXDIM == 2);
      double y[2];
      y[0] = FABS(x[0]);
      y[1] = FABS(x[1]);
      C->D(y, info, cov, v); 
    }
  }
  INFO_TRACE_RETURN("D2", info, cov);
}

void DD_2(double *x, int *info, model *cov, double *v) {
  INFO_TRACE("D2", info, cov);
  assert(everyCoord(isSpaceIsotropic, cov));
  defn *C = DefList + COVNR;// nicht gatternr
  assert(everyCoord(isIsotropic, cov));
  int dim = GATTERXDIM(0); // frueher prevxdim
  if (dim == 1) {
    double y = FABS(*x);

    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht
    C->D2(&y, info, cov, v);// nicht gatternr
  } else {
    //   assert(PREVTOTALXDIM == 2);
    assert(dim == 2);
    system_type *def = DEF;
    if (isIsotropic(def)) {
      double
	x2 = x[0] * x[0],
	t2 = x[1] * x[1],
	r2 = x2 + t2,
	r   = SQRT(r2); 
      
      // (c'(r) * x/r)' = c''(r) * x^2/r^2 + c'(r) [ 1/r - x^2 / r^3]
      C->D2(&r, info, cov, v);// nicht gatternr
      if (r != 0.0) {
	double w;
	C->D(&r, info, cov, &w);
	w /= r;
	*v = (*v - w) * x2 / r2 + w;
      }// else nothing to do
    } else if (equalsSpaceIsotropic(def)) {
      double y[2];
      y[0] = FABS(x[0]);
      y[1] = FABS(x[1]);
      C->D2(y, info, cov, v); // nicht gatternr
    } else BUG;
  }
  INFO_TRACE_RETURN("D2", info, cov);
}

void D3_2(double *x, int *info, model *cov, double *v) {
  INFO_TRACE("D3", info, cov);
  assert(everyCoord(isSpaceIsotropic, cov));
  defn *C = DefList + COVNR;// nicht gatternr
  assert(everyCoord(isIsotropic, cov));
  int dim = GATTERXDIM(0);
  if (dim == 1) {
    double y = FABS(*x);
    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht
    C->D3(&y, info, cov, v);// nicht gatternr
  } else {
    //   assert(PREVTOTALXDIM == 2);
    assert(dim == 2);
    system_type *def = DEF;
    if (isIsotropic(def)) {
      double 
	x2 = x[0] * x[0],
	t2 = x[1] * x[1],
	r2 = x2 + t2,
	r   = SQRT(r2); 
      
      C->D3(&r, info, cov, v);
      if (r != 0.0) {
	double D1, D2,
	  u = x[0] / r,
	  u2 = u * u,
	  w = x[0] / r2;
	C->D(&r, info, cov, &D1);
	C->D2(&r, info, cov, &D2);       
	*v =  *v * u * u2 + 3.0 * (1 - u2) * w * (D2 - D1 / r);
      } 
    } else if (equalsSpaceIsotropic(def)) {
      double y[2];
      y[0] = FABS(x[0]);
      y[1] = FABS(x[1]);
      C->D3(y, info, cov, v); // nicht gatternr
    } else BUG;
  }
  INFO_TRACE_RETURN("D3", info, cov);  
}


void D4_2(double *x, int *info, model *cov, double *v) {
  INFO_TRACE("D4", info, cov);
 assert(everyCoord(isSpaceIsotropic, cov));
  defn *C = DefList + COVNR;// nicht gatternr
  assert(everyCoord(isIsotropic, cov));
  int dim = GATTERXDIM(0);
  if (dim == 1) {
    double y = FABS(*x);
    // iso2iso ueberpruefen !!!
    // dollar + scale // aniso > 0 und dim = 1
    // nach tbm eigentlich auch nicht
    C->D4(&y, info, cov, v);// nicht gatternr
  } else {
    //   assert(PREVTOTALXDIM == 2);
    assert(dim == 2);
    system_type *def = DEF;
    if (isIsotropic(def)) {
      double 
	x2 = x[0] * x[0],
	t2 = x[1] * x[1],
	r2 = x2 + t2,
	r   = SQRT(r2); 
      
      C->D4(&r, info, cov, v);
      if (r != 0.0) {
	double D1, D2, D3,
	  u = x[0] / r,
	  u2 = u * u,
	  u4 = u2 * u2,
	  w = x2 / r;
	C->D(&r, info, cov, &D1);
	C->D2(&r, info, cov, &D2);       
	C->D3(&r, info, cov, &D3);       
	*v = *v * u4 + 6.0 * D3 * u * w * (1 - u2) 
	  + 3.0 * (1 - 6 * u2 + 5 * u4) * (D2  - D1 / r) / r2;
      } 
    } else if (equalsSpaceIsotropic(def)) {
      double y[2];
      y[0] = FABS(x[0]);
      y[1] = FABS(x[1]);
      C->D4(y, info, cov, v); // nicht gatternr
    } else BUG;
  }
  INFO_TRACE_RETURN("D4", info, cov);

}


void inverse2(double *x, model *cov, double *v) {
  defn *C = DefList + COVNR;// nicht gatternr
  C->inverse(x, cov, v);//  nicht gatternr
}
   

void inverse_nonstat2(double *v, model *cov, double *x, double *y){
  defn *C = DefList + COVNR;// nicht gatternr

  C->inverse_nonstat(v, cov, x, y);//  nicht gatternr
}

void inverse_log_nonstat2(double *v, model *cov, double *x, double *y){
  defn *C = DefList + COVNR;// nicht gatternr

  C->nonstat_loginverse(v, cov, x, y);//  nicht gatternr
}
   
   

int struct2(model *cov, model **newmodel) {
  int err;
  errorloc_type errloc_save;
  char *error_location = cov->base->error_location;  

  if (!cov->checked) {
    BUG;
  }
  STRCPY(errloc_save, error_location);
  SPRINTF(error_location, "setting up %.50s", NICK(cov));
 
  err = DefList[COVNR].Struct(cov, newmodel);
   if (newmodel != NULL && (*newmodel) != NULL) {
    SET_CALLING(*newmodel, cov->calling != NULL ? cov->calling : cov);
  }
 
  if (err == NOERROR) STRCPY(error_location, errloc_save);

  RETURN_ERR(err);
}

int init2(model *cov, gen_storage *s){ // s wird durchgereicht!  
  defn *C = DefList + COVNR; //  nicht gatternr
  model
    *calling = cov->calling == NULL ? cov : cov->calling;
  int i,
    err = NOERROR,
    kappas = DefList[COVNR].kappas;
  errorloc_type errloc_save;
  char *error_location = cov->base->error_location;  
  STRCPY(errloc_save, error_location);
   
  //  PrInL++;

  for (i=0; i<kappas; i++) {
    model *param  = cov->kappasub[i];
    if (param != NULL) {
      //      PMI0(param);
      if (isnowRandom(param)) {
	if ((err = INIT_RANDOM(param, 0, s, P(i))) != NOERROR) RETURN_ERR(err);
      } else if (!isnowShape(param) && (err = INIT(param, 0, s)) != NOERROR)
	RETURN_ERR(err);
    }
  }

  if (cov->method == Forbidden) { cov->method = calling->method; }
  
  SPRINTF(error_location, "Initializing %.50s", NICK(cov));
  if (!equalsBernoulliProcess(cov)) {
    switch(cov->frame) {
    case GaussMethodType :
      if (cov->method==SpectralTBM) {
	int nr = COVNR;
	if (cov->calling == NULL && nr != SPECTRAL_PROC_USER &&
	    nr != SPECTRAL_PROC_INTERN) SERR("unexpected value in init2");
      }
      break;
    case InterfaceType: break;
    case BrMethodType: case SmithType: case SchlatherType :
    case PoissonType: case PoissonGaussType: case RandomType:
      cov->origrf = false;
      assert((cov->mpp.moments < 0) xor (cov->mpp.mM != NULL));
      break;
    case TrendType: case EvaluationType: case LikelihoodType: break;
    case NormedProcessType: break;
    default :
      // PMI0(cov);      crash();
      BUG; // ILLEGAL_FRAME;
    }
  }

  if (!cov->initialised && (err = C->Init(cov, s)) != NOERROR)
    goto ErrorHandling;
  
  calling->fieldreturn = cov->fieldreturn;
  
 ErrorHandling :
  //  PrInL--;
  if (err == NOERROR) STRCPY(error_location, errloc_save);
  cov->initialised = err == NOERROR;
  SPRINTF(error_location, "'%.50s'", NICK(calling));//  nicht gatternr   
  RETURN_ERR(err);
}

void do2(model *cov, gen_storage *s){
  //   model *calling = cov->calling == NULL ? cov : cov->calling;
  //
  //  int i,
  //    kappas = DefList[COVNR].kappas;

  // statt nachfolgende Zeilen: siehe init2
  //  for (i=0; i<kappas; i++) {
  //    model *param  = cov->kappasub[i];
  //    if (param != NULL && isnowRandom(param)) DORANDOM(param, P(i));
  //  }

  // printf("do2 hier %s\n", NAME(cov));

  DefList[COVNR].Do(cov, s); // ok

  //  printf("do2 ende hier\n");
  // assert(false);
}


void dorandom2(model *cov, double *v){
 DefList[COVNR].DoRandom(cov, v); // ok
}



void nonstat2statcov(double *x, double *y, int*info, model *cov, double *v) {
  assert(equalsXonly(OWNDOM(0)));
  int dim = PREVTOTALXDIM;
  TALLOC_GATTER(z, PREVTOTALXDIM);
  for (int i=0; i<dim; i++) z[i] = x[i] - y[i];  
  COV(z, info, cov, v);
  END_TALLOC_z;
}


void lognonstat2statcov(double *x, double *y, int*info, model *cov, double *v,
			double *sign) {
  assert(equalsXonly(OWNDOM(0)));
  int dim = PREVTOTALXDIM;
  TALLOC_GATTER(z, PREVTOTALXDIM);
  for (int i=0; i<dim; i++) z[i] = x[i] - y[i];
  LOGCOV(z, info, cov, v, sign);
  END_TALLOC_z;
}


/*
void checkfor(model *cov) {// served by CHECKFOR; nothing to do with CHECK & err
  while (cov!=NULL && !isProcess(cov)) {
    PRINTF("%s ", NAME(cov));
    cov=cov->calling;
  }
  if (cov == NULL) { PRINTF("not process(es)\n"); return; }
  if (OWNISO(0) != PREVISO(0) && OWNISO(0) != ISO_MISMATCH) { APMI(cov); } //
}
*/


sortsofparam SortOf(model *cov,int k, int row, int col, sort_origin origin) {
  defn *C = DefList + COVNR;
  if (C->sortof != NULL) return C->sortof(cov, k, row, col, origin);
  // for MLE
// k non-negative: k-th parameter
// k negative: (i-1)th submodel
  if (k >= C->kappas) BUG;

  return k<0 ? VARPARAM : C->sortof_tab[k];
// k<0: varparam used to indicate that variance is allowed for submodel,
// see recursive call getnaposition; E.g. not allowed for second submodel
// of nsst: the variance parameter is something between scale and variance
}


bool allowedDfalse(model *cov) {
  for(int i=FIRST_DOMAIN; i<= LAST_DOMAINUSER; i++) cov->allowedD[i] = true;
  return false;
}


bool allowedDtrue(model *cov) {
  for(int i=FIRST_DOMAIN; i<=LAST_DOMAINUSER; i++) cov->allowedD[i] = true;
  return true;
}


bool allowedD(model *cov) {
  // assert(!cov->DallowedDone); // 20.10.19 auskommentiert da bei
  // dom=PREVMODEL_D wiederholt aufruft.
  
  // muss zwingend als erstes stehen
  cov->DallowedDone = cov->calling == NULL ? true : cov->calling->DallowedDone;
  //                                  da ParamDepD fehlschlagen kann
  //printf("Dsub: %s\n", NAME(cov));
  assert((int) FIRST_DOMAIN == 0);
  // keine Varianten bei DOMAIN
  defn *C = DefList + COVNR;
  bool *a = cov->allowedD;

  cov->variant=0;
  assert(!isSubModelD(C) || C->Dallowed != NULL);
  if ((C->Dallowed != NULL)) {
    //printf("subi\n");
    model *calling = cov->calling;
    assert(calling != NULL);
    cov->prevloc = LocP(calling);
    return C->Dallowed(cov);
  }

  domain_type dom = DEFDOM(0);
  if (isParamDepD(C) && C->setDI !=NULL && !isFixed(dom) &&
      !C->setDI(cov)){
    cov->DallowedDone = false;
    return  allowedDfalse(cov);
  }
  if (isFixed(dom)) {
    for (int i=(int) FIRST_DOMAIN; i<=(int) LAST_DOMAINUSER; a[i++] = false);
    a[dom] = true;    
    return false;
  }
  return allowedDfalse(cov);
}


bool allowedIfalse(model *cov) {
  //printf("allowedfalse %s\n", NAME(cov));
  bool *I = cov->allowedI;
  for(int i=FIRST_ISOUSER; i<= LAST_ISOUSER; i++) I[i] = true;
  return false;
}


bool allowedItrue(model *cov) {
  //printf("allowedtrue %s\n", NAME(cov));
  bool *I = cov->allowedI;
  for(int i=FIRST_ISOUSER; i<=LAST_ISOUSER; i++) I[i] = true;
  return true;
}


bool allowedI(model *cov) {
  //  printf("all %s\n", NAME(cov));
  // if (cov->IallowedDone) return false; // 20.10.19 auskommentiert wegen wiederholten aufruf bei iso=PREVMODEL_I

  //printf("allowedI:%s %d\n", NAME(cov), cov->zaehler);
  cov->IallowedDone = cov->calling == NULL ? true : cov->calling->IallowedDone;
  //                                  da ParamDepI fehlschlagen kann
  //  printf("Isub: %s (%d)\n", NAME(cov), cov->zaehler);
  assert((int) FIRST_ISOUSER == 0);
  defn *C = DefList + COVNR;
  int variants = C->variants;
  bool *I = cov->allowedI;

  //  printf("hiere %d \n", C->Iallowed != NULL);

 cov->variant=0;
  assert(!isSubModelI(C) || C->Iallowed != NULL);
  if (C->Iallowed != NULL) {
    model *calling = cov->calling;
    assert(calling != NULL);
    cov->prevloc = LocP(calling);
    return C->Iallowed(cov);
  }

  
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  
  isotropy_type iso = DEFISO(0);
  if (isParamDepI(C) && C->setDI!=NULL && !isFixed(iso) && !C->setDI(cov)){
    cov->IallowedDone = false;
    return allowedIfalse(cov);
  }

  if (isFixed(iso)) {
    //    printf("fixed %s\n", ISO_NAMES[iso]);
    I[iso] = true;
    if (equalsUnreduced(iso)) {
      //I[SYMMETRIC] = I[EARTH_SYMMETRIC] = I[SPHERICAL_SYMMETRIC] = false;
      I[CARTESIAN_COORD] =I[EARTH_COORD] = I[SPHERICAL_COORD] = true;
    }
  } else {
    return allowedIfalse(cov);
  }

  for (cov->variant++ ; cov->variant < variants; cov->variant++) {
    // hier nur noch fixe werte zulaessig!
    //   printf("ABBx %d\n", cov->variant);
    iso = DEFISO(0);
    assert(!isParamDepI(C) && !isSubModelI(C));
    assert(isFixed(iso));
    I[iso] = true;  
  }
  
  //  printf("ABB\n");
   cov->variant=0;

  return false;
}


bool allowedIsubs(model *cov, model **sub, int z){
  //  printf("allowed check for %s\n", NAME(cov));
  
  bool *I = cov->allowedI;
  int
    j=0,
    idx_C = FIRST_CARTESIAN,
    idx_E = FIRST_EARTH,
    idx_S = FIRST_SPHERICAL;

  if (z == 0) {
    return allowedItrue(cov);
  }

  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  while (true) {
    if (!allowedI(sub[j])) break;
    if (++j >= z) return allowedItrue(cov);
  }
  MEMCOPY(I, sub[j]->allowedI, sizeof(allowedI_type));
  //  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; i++)  printf("%d \n", I[i]);
  //  printf("%s j=%d\n", NAME(sub[j]), j);  printI(sub[j]);
      
  while (idx_C <= (int) CARTESIAN_COORD && !I[idx_C]) idx_C++;
  while (idx_E <= (int) EARTH_COORD && !I[idx_E]) idx_E++;
  while (idx_S <= (int) SPHERICAL_COORD && !I[idx_S]) idx_S++;
  //PMI(sub[j]);
  //  if (!(idx_C <= (int) CARTESIAN_COORD || idx_E <= (int) EARTH_COORD ||
  //	 idx_S <= (int) SPHERICAL_COORD)) {
    //    PMI(cov);    printI(cov);    printI(cov->sub[0]);
    //    printf("allowIsub %s %d<%d %d<%d %d<%d\n", NAME(cov), idx_C,CARTESIAN_COORD, idx_E, EARTH_COORD ,idx_S,SPHERICAL_COORD);
    //crash();
  //  }
  assert(idx_C <= (int) CARTESIAN_COORD || idx_E <= (int) EARTH_COORD ||
	 idx_S <= (int) SPHERICAL_COORD);
  if (idx_C > CARTESIAN_COORD) {
    if (idx_E > EARTH_COORD) {
      idx_E = idx_S-FIRST_SPHERICAL+FIRST_EARTH;
      for (int i=idx_S; i<=SPHERICAL_COORD; i++)
	I[i-FIRST_SPHERICAL+FIRST_EARTH] = I[i];
    }
  }
    
  for (j++; j<z; j++) {
    assert (sub[j] != NULL);
    if (allowedI(sub[j])) continue;
    bool *subI = sub[j]->allowedI;
    int sub_C = FIRST_CARTESIAN;
    while (sub_C <= (int) CARTESIAN_COORD && !subI[sub_C]) sub_C++;
    int sub_E = FIRST_EARTH;
    while (sub_E <= (int) EARTH_COORD && !subI[sub_E]) sub_E++;
    int sub_S = FIRST_SPHERICAL;
    while (sub_S <= (int) SPHERICAL_COORD && !subI[sub_S]) sub_S++;
 
    if (sub_C <= (int) CARTESIAN_COORD) {
      if (idx_C <= (int) CARTESIAN_COORD) {
	if (equalsVectorIsotropic((isotropy_type) idx_C) && sub_C < idx_C)
	  idx_C = sub_C;
	for (; idx_C < sub_C; I[idx_C++] = false);
	for (int i = idx_C; i<=(int) CARTESIAN_COORD; i++) I[i] |= subI[i];
      }
    } else {
      idx_C = sub_C;
    }
    if (sub_E > EARTH_COORD) {
      if (sub_S <= SPHERICAL_COORD) {
	sub_E = sub_S-FIRST_SPHERICAL+FIRST_EARTH;
	for (int i=sub_S; i<=SPHERICAL_COORD; i++)
	  subI[i-FIRST_SPHERICAL+FIRST_EARTH] = subI[i];
      } else idx_S = sub_S;
    }

    if (sub_S <= (int) SPHERICAL_COORD) {
      if (idx_S <= (int) SPHERICAL_COORD) {
	for (; idx_S < sub_S; I[idx_S++] = false);
	for (int i = idx_S; i<=(int) SPHERICAL_COORD; i++) I[i] |= subI[i];
      } // else {
      //for (; sub_S <=(int) SPHERICAL_COORD; sub_S++) I[sub_S] = subI[sub_S];
      //}
    }
      
    if (sub_E <= (int) EARTH_COORD) {
      if (idx_E <= (int) EARTH_COORD) {
	for (; idx_E < sub_E; I[idx_E++] = false);
	for (int i = idx_E; i<=(int) EARTH_COORD; i++) I[i] |= subI[i];
      } //else {
      //for ( ; sub_E <=(int) EARTH_COORD; sub_E++) I[sub_E] = subI[sub_E];
      // }
    }
    if (idx_C >= (int) CARTESIAN_COORD && idx_E >= (int) EARTH_COORD &&
	idx_S >= (int) SPHERICAL_COORD) break;
  }
  if (idx_C <= CARTESIAN_COORD){
    if (equalsSymmetric((isotropy_type) idx_C)) I[EARTH_SYMMETRIC] = true;
    else I[EARTH_COORD] = I[EARTH_SYMMETRIC] = true;
  }
  if (idx_S <= SPHERICAL_COORD) {
    if (equalsSphericalCoord((isotropy_type) idx_S) &&
	isSymmetric((isotropy_type) idx_C)) // !! not equalsSymmetric!!
      idx_S = SPHERICAL_SYMMETRIC;    
    for (int i=idx_S; i<=SPHERICAL_COORD; i++)
      I[i-FIRST_SPHERICAL+FIRST_EARTH] = I[i] = true;
  } else if (idx_C > CARTESIAN_COORD) {
    if (idx_E >= EARTH_COORD) idx_E = EARTH_SYMMETRIC;
    for (int i=idx_E; i<=EARTH_COORD; i++) I[i] = true;
  }

  //  printf("%s:", NAME(cov));
  //  for (int i=0; i<= EARTH_COORD; i++) {printf("%d.", I[i]);} printf("\n");
  
  return false;
}

