/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2004 Yindeng Jiang & Martin Schlather
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
along with this program; if not, write to the Fre Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/* 

  PRINTING LEVELS
  ===============

  error messages: 1
  forcation : 2
  minor tracing information : 3--5
  large debugging information: >10

*/

/*
  calculate the transformation of the points only once and store the result in
  a register (indexed by ncov) of key, not of s->. Advantages: points have to
  calculated only once for several tries; if zonal anisotropy then specific 
  methods may be called;  
*/

#include "questions.h"
#include "Coordinate_systems.h"
#include "variogramAndCo.h"
#include "Processes.h"
#include "rf_interfaces.h"

#ifdef DO_PARALLEL
#include <omp.h>
#endif

/* 
 in CheckCovariance and other the following dimensions are used:
 xdimOZ       : dimension of the points given by the user.
              value is <= spacedim since in case of isotropy only
	      the distances might be given by the user
 timespacedim : (= kc->dim = key->timespacedim)
              the true dimension for the location (if necessary by 
	      explicite parameter, e.g. in CovarianceFct.)
*/


#define DENS_LOG 0
#define DENS_SEED 1
#define DENS_ENV 2
/*
void density(double VARIABLE_IS_NOT_USED *value, model *cov, double *v) {
KEY_type *K T  = co v->b ase;
  assert(!P0INT(DENS_LOG));
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  char errorloc_save[nErr orLoc];
  int ni = 0,
    err = NOERROR;
  double *res;
  locati on_type *l oc = LocPrev(cov);
  Long vdimtot = (Long) loc->totalpoints * VDIM0;
  assert(VDIM0 == VDIM1);

  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array

  STRCPY(errorloc_save, K T->err or_location);

  PutRN Gstate();
  ERR("stop : ni nae Zei falsch");
  double simu_seed = GL OBAL_UTILS->basic.seed + (ni - 1);
  addVariable((char*) "seed", &simu_seed, 1, 1, PENV(DENS_ENV)->sexp);
  eval(PLANG(DENS_SEED)->sexp, PENV(DENS_ENV)->sexp);
  GetRN Gstate();

  SPRINTF(K T->erro r_location, "%.50s %d", errorloc_save, ni);
 
  assert(cov->Sgen != NULL);
  
  NotProgrammedYet("density");

  //  if (COVNR == DENSITY) LOGDENSITY(sub, cov->Sgen); 
  //  else if (COVNR == PROBAB) PROBAB(sub, cov->Sgen);

  if (sizeof(double) == sizeof(double) && false) {
    MEMCOPY(res, cov->rf, sizeof(double) * vdimtot);
  } else {
    int i; for (i=0; i<vdimtot; i++) {
      res[i] = cov->rf[i];
    }
  }
  
  if (!sub->simu.ac tive) 
    GERR1("could not perform multiple simulations. Is '%.50s == FALSE'?",
	  general[GENERAL_STORING]);

 ErrorHandling: 
 
  PutRN Gstate();
  
  if (err > NOERROR) {
    XERR(err);
  }
}


int check_density(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  lo cation_type *l oc = L ocPrev(cov);
  int j, err, frame;
  isotropy_type iso;
  Types type;
  //  bool vdim_close_together = glob al->general.vdim_close_together;

  ASSERT_LOC_GIVEN;

 
  cov->simu.expected_number_simu = gl obal->general.expected_number_simu;
  if (cov->simu.expected_number_simu > 1 && 
      gl obal->general.expected_number_simu <= 1)
    SERR("expected number of simulations inconsistent");

  glo bal->internal.sto red_init = glob al->general.storing || 
    gl obal->general.expected_number_simu > 1;  
  
  if (cov->key == NULL) {
    domain_type dom = KERNEL;
    if (isProcess(sub)) {
      frame = GaussMethodType;
      type = ProcessType;
      iso = UNREDUCED;
    } else {
      frame = V ariogramType;
      type = PosDefType;
      iso = S YMMETRIC;
    }
    if (cov->frame = = Any Type) frame = Any Type;

    err = ERR ORTY PECONSISTENCY;

    for (j=0; j<=2; j++) {
      if ((Type Consistency(type, sub) && 
	   (err = C HECK(sub, loc->timespacedim, OWNXDIM(0), type, 
	      dom, iso, cov->vdim, frame)) == NOERROR) 
	  || isProcess(sub)) break;

      if (j==0) type = VariogramType;
      else {
	type = TrendType;
	dom = X ONLY;
	iso = CARTESIAN_COORD;
      }
    }

    if (err != NOERROR) RETURN_ERR(err);
  } else {
    BUG;
    frame = fram e_o f_p rocess(SUBNR);
    if (fram e == R OLE_FAILED) BUG;

    if ((err = C HECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType,
		     X ONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD : PREVISO(0),
		     cov->vdim, frame)) != NOERROR) {
      RETURN_ERR(err);
    }
  }

  setbackward(cov, sub);  
  int subvdim = sub->vdim[0];
  VDIM0=subvdim; 
  VDIM1=sub->vdim[1]; 

  if (cov->q == NULL) {
    QALLOC(1);
    cov->q[0] = 1.0; 
  }

  RETURN_NOERROR;
}


int struct_density(model *cov, model VARIABLE_IS_NOT_USED  **newmodel){
  model *next = cov->sub[0],
    *sub = next;
  locati on_type *l oc = Lo cPrev(cov);
  int err,
    subframe = R OLE_FAILED,
    nr = NEXTNR;

  //APMI(next);

  if (isVariogram(next) || is Trend(next)) {
    if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
    addModel(&(cov->key), GAUSSPROC);
    sub = cov->key;
    
    if ((err = C HECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType,
		     X ONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD : PREVISO(0), 
		     cov->vdim, GaussMethodType)) != NOERROR) {
      RETURN_ERR(err);
    }
    subframe = GaussMethodType;    
  } else if (is BernoulliProcess(next)) subframe = FRAME_ BERNOULLI;
  else if (isGaussMethod(next)) subframe = GaussMethodType;   
  else if (isBrMethod(next)) subframe = BrMethodType;    
  else if (nr = = POISSONPROC) subframe = Poisson;
  else if (nr = = SCHLATHERPROC) subframe = SchlatherType;
  else if (nr = = EXTREMALTPROC) subframe = SchlatherType;
  else if (nr = = SMIT HPROC) subframe = SmithType;
  else {
    ILLEGAL_FRAME;
  }
  sub->frame = subframe; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C-HECK (was auf jeden Fall gemacht werden muss)

  cov->simu.act ive = next->simu.act ive = false; 
  sub->simu.expected_number_simu = cov->simu.expected_number_simu;

  if (PL >= PL_DETAILS) PRINTF("Struct Density\n");

  if ((err = STRUCT(sub, NULL)) != NOERROR) { RETURN_ERR(err); }
  if (PL >= PL_DETAILS) PRINTF("Checking Density\n");


  assert(cov->Sgen == NULL);
  NEW_STORAGE(gen);

  if (!sub->initialised) {
    if (PL >= PL_DETAILS) PRINTF("Struct Density C\n");
   
    if (//cov->key ! = NULL && // bei sub= =next waere der falsche frame gesetzt
	// irgenwie komisch, cov->key abzufragen und check(sub aufzurufen
	// aufgrund von Beispiel in RTpoisson geloescht
	(err = C HECK(sub, loc->timespacedim, OWNXDIM(0), ProcessType,
		     PREVDOM(0), PREVISO(0), cov->vdim,
		     subframe)) != NOERROR) {
      //
      RETURN_ERR(err);
    }
    
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nStruct Density (%s, #=%d), after 2nd check:",
	     NICK(sub), sub-> gatter nr); // OK
      PMI(sub); // OK
    }
 
    
    if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
      RETURN_ERR(err); 
    }
  }
 
  cov->initialised = true;
  ReturnOtherField(cov, sub);

  cov->simu.act ive = sub->simu.act ive = true; 
  RETURN_NOERROR;
}

void range_density(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  booleanRange(DENS_LOG);
}
*/

 

int alloc_fctn(model *cov, int xdim, int rowscols) { 
  /*
    
    alloc_fctn:
    direct.cc, gausslikeli(Fct nIntern, CovarianceMatrix), startGetNset (StandardCovMatrix), RMS (covmatrixS)

    che ck_fct_intern called by ch eck_cov_intern, che ck_fctn (RFFct n), check_predict (RFpredict),
    che ck_cov_int called by check_vario (RFvario), ch eck_cov (struct_cov, RFCov, RFPseudoVario)
    che ck_covmatri called by struct_cov, RFCovMatrix
    

    fctn->call (nur variogramAndXCo & rf_interfaces)
 

    as such: gridlen, start, end, nx, x, xstart, inc, endy, startny, ny, ystart, incy, z, Val, cross, C0x, C0y
    

   */
  int dimP1 = xdim + 1; // safety -- some need it

  if (cov->Sfctn != NULL) fctn_DELETE(&(cov->Sfctn));
  NEW_STORAGE(fctn);
  getStorage(fctn ,   fctn); 

  //  printf("dimP1 = %d\n", dimP1);

  fctn->rowscols = rowscols;
  if (  // grid
      (fctn->start = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (fctn->end = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (fctn->nx = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (fctn->startny = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (fctn->endy = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
      (fctn->ny = (int*) CALLOC(dimP1, sizeof(int))) == NULL ||
     
      (fctn->x = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (fctn->xstart = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (fctn->inc = (double*) CALLOC(dimP1, sizeof(double))) == NULL||
      (fctn->y = (double*) CALLOC(dimP1, sizeof(double))) == NULL || // auch fct
      (fctn->ystart = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
      (fctn->incy = (double*) CALLOC(dimP1, sizeof(double))) == NULL ||
   
      // CovarianceMatrix only
      (fctn->cum = (int*) CALLOC(dimP1, sizeof(int))) == NULL
	)
    RETURN_ERR(ERRORMEMORYALLOCATION);

  // fuer CovVario:
  if (rowscols > 0 &&
      ((fctn->cross= (double*) CALLOC(rowscols, sizeof(double))) == NULL ||
       //   (fctn->Val= (double**) CALLOC(rowscols, sizeof(double*))) == NULL ||  
       (fctn->C0x  = (double*) CALLOC(rowscols, sizeof(double))) == NULL || 
       (fctn->C0y  = (double*) CALLOC(rowscols, sizeof(double))) == NULL
       ))
    RETURN_ERR(ERRORMEMORYALLOCATION);
  

  RETURN_NOERROR;
}



#define SIMU_CHECKONLY 0
#define SIMU_SEED 1
#define SIMU_ENV 2

void simulate(double VARIABLE_IS_NOT_USED *X, int *N, model *cov, double *v){
  globalparam *global = &(cov->base->global);
  utilsparam *global_utils = &(cov->base->global_utils);
  assert(!P0INT(SIMU_CHECKONLY));
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  errorloc_type errorloc_save;
  char
    format[20],
    back[]="\b\b\b\b\b\b\b\b\b\b\b", 
    prozent[]="%",
    pch = global->general.pch;
  int ni, digits, //n, 
    nn,
    err = NOERROR,
    each = 0;
  double *res,
    realeach=0.0;
  simu_storage *simu = NULL;  
  finaldofct finalDo = DefList[SUBNR].FinalDo;
  Long vdimtot = (Long) Loctotalpoints(cov) * VDIM0;
  assert(VDIM0 == VDIM1);
  char *error_location = cov->base->error_location;  
 
  // die folgenden Zeilen duefen nicht abgeaendert werden!!
  // Grund; *N wird eventuell durch Koordinatentrafo veraendert
  // so dass bei v=NULL ohne Umweg ueber gatter aufgerufen wird
  // und der korrekte Wert in cov->q gespeichert wird, der dann
  // spaeter wieder ausgelesen wird und nn zugeordnet wird.
  if (v == NULL) {
    cov->q[cov->qlen - 1] = *N;
    return; // EvaluateModel needs information about size
  //                      of result array
  } 
  nn = (int) cov->q[cov->qlen - 1];

  STRCPY(errorloc_save, error_location);
  simu = &(cov->simu);
  if (!simu->active) {
    err=ERRORNOTINITIALIZED; goto ErrorHandling;
  }

  if (nn>1 && pch != '\0') {
    if (pch == '!') {
      digits = (nn<900000000) 
	? 1 + (int) TRUNC(LOG((double) nn) / LOG(10.0)) : 9;
      back[digits] = '\0';
      each = (nn < 100) ? 1 :  nn / 100;
      SPRINTF(format, "%.2ss%.2s%dd", prozent, prozent, digits);
    } else if (pch == '%') {
      back[4] = '\0';
      realeach = (double) nn / 100.0;
      each = (nn < 100) ? 1 : (int) realeach;
      SPRINTF(format, "%.2ss%.2s%dd%.2ss", prozent, prozent, 3, prozent);
    } else each = 1;
  } else each = nn + 1;
  res = v;

  sub->simu.pair = false;
  for (ni=1; ni<=nn; ni++, res += vdimtot) {
#ifdef DO_PARALLEL
    //omp_set_num_threads(1);
#endif

    if (global_utils->basic.seed != NA_INTEGER &&
	(nn > 1 || global->general.seed_incr != 0)) {
      assert(!PisNULL(SIMU_SEED) && !PisNULL(SIMU_ENV));
      PutRNGstate(); // PutRNGstate muss vor UNPROTECT stehen!! -- hier irrelevant
      // if seed is increased by (ni -1) only, a row of simulations
      // of max-stable fields will look very much the same !!
      int simu_seed = global_utils->basic.seed + (ni - 1) *
	global->general.seed_sub_incr
	+ nn * global->general.seed_incr;

      // printf("seed=%d %d %d\n",simu_seed, GL OBAL_UTILS->basic.seed,global->general.seed_incr);
      
      addIntVariable((char*) "seed", &simu_seed, 1, 1, PENV(SIMU_ENV)->sexp);
      eval(PLANG(SIMU_SEED)->sexp, PENV(SIMU_ENV)->sexp);
      GetRNGstate(); // start interner modus
    }

    
    SPRINTF(error_location, "%.50s %d", errorloc_save, ni); 
    
    if (PL > 0) {
      if (ni % each == 0) {
	if (pch == '!')  
	  PRINTF(format, back, ni / each);
	else if (pch == '%')
	  PRINTF(format, back, (int) (ni / realeach), prozent);
      else PRINTF("%c", pch);
      }
    }
  
    assert(cov->Sgen != NULL);
   
    DO(sub, cov->Sgen);

    //    printf("ni=%d vdimtot=%d %ld sub=%ld\n", ni, vdimtot, res, sub->rf);
    //if(false) for (int k=0; k<vdimtot; k++)  if (sub->rf[k] <= 0.0) APMI(sub);

    
#ifdef DO_PARALLEL
    if (global_utils->basic.cores>1)
      omp_set_num_threads(global_utils->basic.cores);
    R_CheckUserInterrupt();
    if (global_utils->basic.cores>1) omp_set_num_threads(1); // !! wichtig,
    // auch fuer nachfolgende DO()!!
#endif    
    
    MEMCOPY(res, cov->rf, sizeof(double) * vdimtot);
    // for (int i=0; i<vdimtot; i++) res[i] = cov->rf[i];
  
    if (!sub->simu.active) 
      GERR1("could not perform multiple simulations. Is '%.50s == FALSE'?",
	    general[GENERAL_STORING]);
  }
  
  if (finalDo != NULL) {
    finalDo(sub, v, nn, cov->Sgen);
  }

#ifdef DO_PARALLEL
  omp_set_num_threads(global_utils->basic.cores);
#endif    

  if (nn>1 && pch != '\0') {
    if (pch == '!' || pch == '%') PRINTF("%s", back);
    else PRINTF("\n");
  }

  assert(simu != NULL);
  sub->simu.active = simu->active = sub->simu.active && global->general.storing;

 ErrorHandling:

  //  PMI(cov);
   
  if (err > NOERROR) {
    if (simu != NULL) sub->simu.active = simu->active = false;
    XERR(err);
  }
}


int check_simulate(model *cov) {
  globalparam *global = &(cov->base->global);
  utilsparam *global_utils = &(cov->base->global_utils);
 model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  int j, d, 
    err = NOERROR,
    dim = Loctsdim(cov);
  isotropy_type iso;
  Types type, frame;
  bool vdim_close_together = global->general.vdim_close_together;
  
  cov->initialised = false;
  ASSERT_LOC_GIVEN;
  kdefault(cov, SIMU_CHECKONLY, false);

  cov->simu.expected_number_simu = global->general.expected_number_simu;
  if (cov->simu.expected_number_simu > 1 && 
      global->general.expected_number_simu <= 1)
    SERR("expected number of simulations inconsistent");
 
  if (global->general.seed_incr != 0 && global_utils->basic.seed == NA_INTEGER){
    SERR("'seed' must be set if 'seed_incr' is different from 0.");
  }

  cov->base->stored_init = global->general.storing || 
    global->general.expected_number_simu > 1;  
  
  if (cov->key == NULL) {
    domain_type dom;

    //    PMI(sub);
    //    printf("%d\n", isProcess(sub));
    
    if (isProcess(sub)) {
      dom = XONLY; // 5/2015; voher beides Kerne
      frame = InterfaceType, //InterfaceType; // ProcessType; // Interface
      type = ProcessType;
      // iso = CoordinateSystemOf(PREVISO(0));
    } else {
      dom = isAnyIsotropic(OWNISO(0)) ? XONLY : KERNEL;
      frame = EvaluationType;
      type = PosDefType;
    }
    iso = PREVISO(0);
    assert(equalsCoordinateSystem(iso) || isAnyIsotropic(iso));

    if (hasAnyEvaluationFrame(cov)) BUG;
    // if (cov->frame = = Any Type) frame = Any Type;

    int errold = CERRORTYPECONSISTENCY;
    for (j=0; j<=2; j++) {
       //      PMI0(cov);
      // printf("\nj=%d ownxdim = %d %s\n", j, OWNXDIM(0), TYPE_NAMES[type]);
      // BUG;
     if ( /// auskommenieren ? 8.8.17 ?? Type Consistency(type, sub, 0) &&
	 (err = CHECK(sub, dim, OWNXDIM(0), type, 
		       dom, iso, cov->vdim, frame)) == NOERROR) {
	break;
      }
     
      if (isProcess(sub)) {// muss als zweites stehen !!
	break;
      }

      if (j==0) {
	errold = err; // Fehler(meldungen) werden abstrus, falsch model falsch
	type = VariogramType;
	errold = err;
      } else {
	type = TrendType;
	dom = XONLY;
	iso = PREVISO(0);
      }
    }
    if (err != NOERROR) {
      err = j==0 ? err : errold;
      RETURN_ERR(err);
    }
  } else {
    frame = InterfaceType;
    assert(equalsCoordinateSystem(PREVISO(0)));

    if ((err = CHECK(sub, dim, OWNXDIM(0), ProcessType, XONLY, PREVISO(0),
		     cov->vdim, frame)) != NOERROR) {
      RETURN_ERR(err);
    }
  }
  
  setbackward(cov, sub);  
  int subvdim = sub->vdim[0];
  VDIM0=subvdim; 
  VDIM1=sub->vdim[1]; 

  if (cov->q == NULL) {
    int len=1,
      tsdim = Loctsdim(cov);
    bool grid = Locgrid(cov);
    coord_type gr = Locxgr(cov);
    if (grid) len += tsdim; else len++;
    if (subvdim > 1) len++;
    QALLOC(len);
    cov->q[--len] = 1.0; // number of simulations
    if (subvdim > 1 && !vdim_close_together) cov->q[--len]=(double) subvdim;
    if (grid) {
      for (d=tsdim - 1; d>=0; d--) 
	cov->q[--len] = gr[d][XLENGTH];
    } else {
      cov->q[--len] = Loctotalpoints(cov);
    }
    if (subvdim > 1 && vdim_close_together) {
      assert(len == 1);
      cov->q[--len]=(double) subvdim;
    }
    assert(len==0);
  }

  cov->initialised = true;
  RETURN_NOERROR;
}



int struct_simulate(model *cov, model VARIABLE_IS_NOT_USED  **newmodel){  
  model *next = cov->sub[0],
    *sub = next;
  int err;
  Types subframe = InterfaceType;

  if (isnowVariogram(next) || isnowTrend(next)) {
    if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
    assert(cov->key->calling == next->calling); 
    addModelKey(cov, isnowVariogram(next) ? GAUSSPROC : SHAPE_FCT_PROC);
    sub = cov->key;

    if ((err = CHECK(sub, Loctsdim(cov), OWNXDIM(0), ProcessType, XONLY, 
		     isCartesian(PREVISO(0)) ? CARTESIAN_COORD :PREVISO(0), //??
		     cov->vdim, InterfaceType)) != NOERROR) {
       RETURN_ERR(err);
    }    
  }
  
  sub->frame = subframe; // ansonsten muesste hier C-HECK aufgerufen werden
  // hoffentlich gut gehende Abkuerzung, dass S-TRUCT aufgerufen wird,
  // und danach C-HECK (was auf jeden Fall gemacht werden muss)

  cov->simu.active = next->simu.active = false; 
  sub->simu.expected_number_simu = cov->simu.expected_number_simu;
  if (P0INT(SIMU_CHECKONLY)) RETURN_NOERROR;

  if (PL >= PL_DETAILS) { PRINTF("Struct Simulate\n"); }
 
  if ((err = STRUCT(sub, NULL)) != NOERROR) {
    // printf(">>> %s\n", cov->base->error_location);
    RETURN_ERR(err);
  }
   
  if (PL >= PL_DETAILS) { PRINTF("Checking Simulate\n");}

  NEW_STORAGE(gen);    
  if (!sub->initialised) {
    if (PL >= PL_DETAILS) { PRINTF("Struct Simulate C\n");}
   
    if (//cov->key ! = NULL && // bei sub= =next waere der falsche frame gesetzt
	// irgenwie komisch, cov->key abzufragen und check(sub aufzurufen
	// aufgrund von Beispiel in RTpoisson geloescht
	(err = CHECK_PASSTF(sub, ProcessType, cov->vdim[0], subframe))
	!= NOERROR)
      //	(err = CHECK(sub, loc->timespacedim, cov->xdimown, ProcessType,
      //	     cov->domprev, cov->isoprev, cov->vdim,
      //	     subframe)) != NOERROR)
      RETURN_ERR(err);   
    
    if (PL >= PL_DETAILS) {
      PRINTF("\n\nStruct Simulate (%s, #=%d), after 2nd check:",
	     NICK(sub), SUBNR); // OK
      PMI(sub); // OK
    }
 
    if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
      RETURN_ERR(err); 
    }
  }
 
  cov->initialised = true;
  ReturnOtherField(cov, sub);
  
  //  PMI(sub);
  assert( sub->simu.active );
  cov->simu.active = sub->simu.active; 

  RETURN_NOERROR;
}

void range_simulate(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  booleanRange(SIMU_CHECKONLY);

#define simu_n 2
  int idx[simu_n] = {SIMU_SEED, SIMU_ENV};
  for (int i=0; i<simu_n; i++) {
    int j = idx[i];
    range->min[j] = RF_NAN;
    range->max[j] = RF_NAN;
    range->pmin[j] = RF_NAN;
    range->pmax[j] = RF_NAN;
    range->openmin[j] = false;
    range->openmax[j] = false; 
  }
}

//void do_simulate(model *cov, gen_storage *s){
//  assert(false);
//}

/////////////////////////////////////////




double GetPriors(model *cov) {
  defn *C = DefList + COVNR;
  model *kap;
  int 
    kappas = C->kappas,
    nsub = cov->nsub;
  double v,
    logli = 0.0;
  DEFAULT_INFO(info);
  for (int i=0; i<kappas; i++) {
    if ((kap = cov->kappasub[i]) != NULL) {
      if (isnowRandom(kap)) {
	if (C->kappatype[i] < LISTOF) {
	  // PMI(kap);
	  assert(C->kappatype[i] == REALSXP);
	  VTLG_DLOG(P(i), info, kap, &v);
	} else if (C->kappatype[i] < LISTOF + LISTOF) { // not tested yet
	  NotProgrammedYet("hierachical models for multiple data sets");
	  assert(C->kappatype[i] == LISTOF + REALSXP);
	  int
	    store = cov->base->set,
	    end = cov->nrow[i];
	  for (int set = 0; set<end; set++) {
	    cov->base->set = set;
	    VTLG_DLOG(LP(i), info, kap, &v);
	    logli += v;
	  }
	  cov->base->set = store;
	} else BUG;
	logli += v;
      }
      logli += GetPriors(kap);
    }
  }
  
  for (int i=0; i<nsub; i++) {
    logli += GetPriors(cov->sub[i]);
  }

  return logli; 
}





/* *************** */
/* LOG LIKELIHOOD */
/* *************** */

 void kappalikelihood(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
   *nr = *nc = i <= LIKELIHOOD_DATA ? SIZE_NOT_DETERMINED
    : i <= LIKELIHOOD_LAST ? 1 : -1;
}


void likelihood(double VARIABLE_IS_NOT_USED *x, int *info,
		model *cov, double *v) {
  model *process = cov->key == NULL ? cov->sub[0] : cov->key;
  if (v == NULL) {
    // likelihood_storage *L = process->Slikelihood;
    GETSTORAGE(L ,  process,   likelihood); 
    assert(L != NULL);
    int //betas = L->cum_n_betas[L->fixedtrends],
      vdim = process->vdim[0];
    likelihood_facts *facts = &(L->facts);
    listoftype *datasets = L->datasets;
    int
      store = cov->base->set,
      betatot = L->cum_n_betas[L->fixedtrends],
      all_betatot = betatot;
    cov->base->set = 0;
    if (L->betas_separate)
      all_betatot *= NCOL_OUT_OF(datasets) / vdim;
    cov->q[0] = 1 + facts->globalvariance + all_betatot;
    cov->q[1] = 1;
    cov->base->set = store;
    return; 
  }

  assert(isProcess(process));
  VTLG_DLOG(NULL, info, process, v); 

  assert(process->key == NULL);
  if (isProcess(process->sub[0]))//besser: EmptySubOK(MODELNR(process->sub[0])))
    ERR("processes as priors not allowed yet");
  *v += GetPriors(process->sub[0]); 
 }

 int check_likelihood_linear(model *cov) {
  assert(cov->prevloc != NULL);
  bool distances = LocDist(cov); // ehem Prev
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  int err = CERRORTYPECONSISTENCY,
    logdim = Loctsdim(cov),
    dim = distances ? LocxdimOZ(cov) : logdim;
  domain_type dom = distances == 1 && dim == 1 ? XONLY : KERNEL;
  Types frame = LikelihoodType;
    // SUBNR == PLUS || isProcess(sub) ? LikelihoodType : EvaluationType;

  cov->base->rawConcerns = COVNR == PREDICT_CALL ? noConcerns
    : onlyOneDataSetAllowed;

  
  if (isProcess(sub)) {
    err = CHECK(sub, logdim, dim, ProcessType, dom,
		distances ? ISOTROPIC : UNREDUCED,
		cov->vdim, frame);
  } else {
    err = CHECK(sub, logdim, dim, PosDefType, dom,
		distances ? ISOTROPIC : CoordinateSystemOf(PREVISO(0)),
		cov->vdim,
		frame);
    //    OnErrorStop(err, sub);
    //    APMI0(sub);
    
    
    if (err != NOERROR) {
      err = CHECK(sub, logdim, dim, NegDefType, dom,
		  distances ? ISOTROPIC : SymmetricOf(PREVISO(0)),
		  cov->vdim,
		  frame);
    }
  }
  if (err != NOERROR) RETURN_ERR(err);

  setbackward(cov, sub);  
  VDIM0=sub->vdim[0]; 
  VDIM1=sub->vdim[1]; 

  if (cov->q == NULL) {
    QALLOC(2);
    cov->q[0] = LoctotalpointsY(cov);
  }
  cov->q[1] = VDIM0;

  RETURN_ERR(err);
}


int check_likelihood(model *cov) {
  globalparam *global = &(cov->base->global);
  int err,
    sets = LocSets(cov);
  if ((err = check_likelihood_linear(cov))) RETURN_ERR(err);
  
  kdefault(cov, LIKELIHOOD_NA_VAR, global->fit.estimate_variance);
  kdefault(cov, LIKELIHOOD_BETASSEPARATE, false);
  if (P0INT(LIKELIHOOD_BETASSEPARATE)) BUG; // to estimate beta 
  // for each independent repetition within a data set is a feature
  // that has been implemented, but currently it is not used
  kdefault(cov, LIKELIHOOD_IGNORETREND, false);
  if (PisNULL(LIKELIHOOD_DATA)) BUG;
  for (int set = 0; set<sets; set++) {
    cov->base->set = set;
    Long 
      datatot = LNROW(LIKELIHOOD_DATA) * LNCOL(LIKELIHOOD_DATA),
      totpts = LoctotalpointsY(cov),
      vdimtot = totpts * VDIM0;
    int
      repet = datatot / vdimtot;

    //    PMI0(cov);
    //    printf("repet = %d %d %d %d\n", repet,  vdimtot , VDIM0,datatot);
    
    if (repet * vdimtot != datatot || repet == 0)  {
      GERR("data and coordinates do not match");
    }
    LNRC_(LIKELIHOOD_DATA, nrow) = totpts;
    LNRC_(LIKELIHOOD_DATA, ncol) = datatot / totpts;
  }

 ErrorHandling:
  RETURN_ERR(err);
}

int struct_likelihood(model *cov, 
		      model VARIABLE_IS_NOT_USED  **newmodel){
  //return struct_rftrend(cov, newmodel);
  assert(cov->key == NULL);
  model *sub = cov->sub[0];
  int err = NOERROR;

  if (isnowVariogram(sub)) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) RETURN_ERR(err);
    addModelKey(cov, GAUSSPROC); // struct_gauss_logli calls alloc_fctn
    //                              so that Sfctn is allocated at GAUSSPROC
    sub = cov->key;
    if ((err = CHECK(sub, Loctsdim(cov), OWNXDIM(0), ProcessType, XONLY, 
		     isCartesian(PREVISO(0)) ?
		     (LocDist(cov)  && LocxdimOZ(cov) == 1
		      ? ISOTROPIC : CARTESIAN_COORD)
		     : PREVISO(0),
		     cov->vdim, LikelihoodType)) != NOERROR) {
      RETURN_ERR(err);
    }
  } else sub->frame = LikelihoodType;

  if (!isnowProcess(sub)) 
    SERR1("'%.50s' can be calculated only for processes.", NICK(cov));
  
  if ((err = STRUCT(sub, NULL)) != NOERROR) RETURN_ERR(err);
  assert(sub->Sfctn != NULL);
  NEW_STORAGE(gen);
  if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}


void range_likelihood(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[LIKELIHOOD_DATA] = RF_NEGINF;
  range->max[LIKELIHOOD_DATA] = RF_INF;
  range->pmin[LIKELIHOOD_DATA] = - 1e8;
  range->pmax[LIKELIHOOD_DATA] = 1e8;
  range->openmin[LIKELIHOOD_DATA] = true;
  range->openmax[LIKELIHOOD_DATA] = true;
  
  booleanRange(LIKELIHOOD_NA_VAR);
  booleanRange(LIKELIHOOD_BETASSEPARATE);
  booleanRange(LIKELIHOOD_IGNORETREND);
}



void linearpart(double VARIABLE_IS_NOT_USED *x, INFO, model  VARIABLE_IS_NOT_USED *cov, double  VARIABLE_IS_NOT_USED *v) { 
  BUG; // darf nicht aufgerufen werden
}

SEXP get_linearpart(SEXP model_reg, SEXP Set){
   int cR = INTEGER(model_reg)[0];
   set_currentRegister(cR);
  if (cR < 0 || cR > MODEL_MAX) BUG;	
  model *cov = KEY()[cR];			
  cov = cov->key != NULL ? cov->key : cov->sub[0];
  //  TREE(cov);
  if (COVNR == GAUSSPROC) 
    return gauss_linearpart(model_reg, Set);
  else BUG;
  //  TREE(cov); 
}


int check_linearpart(model *cov) {
  if (LocDist(cov)) 
    SERR1("'%.50s' may not be used when distances are given", NAME(cov));
  return check_likelihood_linear(cov);
}
  



int struct_linearpart(model *cov, 
    model VARIABLE_IS_NOT_USED  **newmodel){
  assert(cov->key == NULL);
  model *sub = cov->sub[0];
  int err = NOERROR;
 
  if (isnowVariogram(sub)) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) RETURN_ERR(err);
    addModelKey(cov, GAUSSPROC);
    sub = cov->key;
    
    if ((err = CHECK(sub, Loctsdim(cov), OWNXDIM(0), ProcessType, XONLY, 
		     isCartesian(PREVISO(0)) ?
		     (LocDist(cov) && LocxdimOZ(cov) == 1 ? ISOTROPIC
		      :CARTESIAN_COORD) : PREVISO(0),
		     cov->vdim, LikelihoodType)) != NOERROR) {
      RETURN_ERR(err);
    }
  } else sub->frame = LikelihoodType;
  
  if (!isnowProcess(sub)) 
    SERR1("'%.50s' can be calculated only for processes.", NICK(cov));


  if ((err = STRUCT(sub, NULL)) != NOERROR) RETURN_ERR(err);
  // likelihood_storage *L = sub->Slikelihood;
  GETSTORAGE(L , sub,  likelihood); 
  if (L == NULL) RETURN_ERR(ERRORFAILED);
  if (L->dettrend_has_nas || L->fixedtrend_has_nas) {
    warning("NAs detected in '%.20s'; hence zero's introduced", NICK(cov)); 
  }
  RETURN_ERR(err);
}




#define EVALDISTR_X 0 // d
#define EVALDISTR_Q 1 // p
#define EVALDISTR_P 2 // q
#define EVALDISTR_N 3 // r
#define EVALDISTR_DIM 4 // r

void kappa_EvalDistr(int i, model VARIABLE_IS_NOT_USED *cov, 
		     int *nr, int *nc){
  *nc = *nr = i <= EVALDISTR_N ? SIZE_NOT_DETERMINED
    : i == EVALDISTR_DIM ? 1 : -1;
}

void EvalDistr(double VARIABLE_IS_NOT_USED *X, int *info, model *cov,
	       double *v){
  errorloc_type errorloc_save;
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  double  *xqp;
  int i, j,
    dim = ANYDIM,
    n = (int) (cov->q[cov->qlen - 1]);

  if (v==NULL) return; // EvaluateModel needs information about size
  STRCPY(errorloc_save, cov->base->error_location);

  if (!PisNULL(EVALDISTR_X)) { // d
    xqp = P(EVALDISTR_X);
    for (j=i=0; i<n; i++, j+=dim) VTLG_D(xqp + j, info, sub, v + i);
  } else if (!PisNULL(EVALDISTR_Q)) { // p
    xqp = P(EVALDISTR_Q);
    for (j=i=0; i<n; i++, j+=dim) VTLG_P(xqp + i, info, sub, v + j);
  } else if (!PisNULL(EVALDISTR_P)) {// q
    xqp = P(EVALDISTR_P);
    for (j=i=0; i<n; i++, j+=dim) VTLG_Q(xqp + j, sub, v + i);
  } else if (!PisNULL(EVALDISTR_N)) {// r
    xqp = P(EVALDISTR_N);
    for (j=i=0; i<n; i++, j+=dim) {
      VTLG_R(NULL, info, sub, v + j);    
    }
  } else BUG;
}


int check_EvalDistr(model *cov) {
  defn *C = DefList + COVNR;
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  int size_idx,  err, //type,  nr = SUBNR,
    *nrow = cov->nrow,
    *ncol = cov->ncol,
    dim = ANYDIM, 
    zaehler=0;

  
  if (cov->q == NULL) {
    int nn = 1; // !! 1 obwohl 2 reserviert werden !! Wg EvaluateModel
    if (dim > 1 && 
	((!PisNULL(EVALDISTR_N) && P0(EVALDISTR_N) > 1) ||
	 (!PisNULL(EVALDISTR_Q) && P0(EVALDISTR_Q) > 1)))
      nn++;
    QALLOC(nn + 1); 
    cov->qlen = nn;
   /*  cov->qlen = 1; // !! 1 obwohl 2 reserviert werden !! Wg EvaluateModel
    if (dim > 1 && 
	((!PisNULL(EVALDISTR_N) && P0(EVALDISTR_N) > 1) ||
	 (!PisNULL(EVALDISTR_Q) && P0(EVALDISTR_Q) > 1)))
      cov->qlen++;
      cov->q = (double *) MALLOC(sizeof(double) * (cov->qlen + 1)); // QALLOC NOT APPROPRIATE
   */

    cov->q[0] = dim; // is overwritten if not a matrix is returned
    size_idx = cov->qlen - 1;

    if (PisNULL(EVALDISTR_N)) {
      if (!PisNULL(EVALDISTR_X)) { // d
	if (dim > 1 && nrow[EVALDISTR_X] != dim) 
	  SERR2("dimenson of '%.50s' does not match '%.50s' ", 
		C->kappanames[EVALDISTR_X], C->kappanames[EVALDISTR_DIM]);
	cov->q[size_idx] = nrow[EVALDISTR_X] * ncol[EVALDISTR_X] / dim;
	zaehler++;
      }
      if (!PisNULL(EVALDISTR_Q)) { // p
	if (dim > 1 && nrow[EVALDISTR_Q] != dim) 
	  SERR2("dimension of '%.50s' does not match '%.50s' ", 
		C->kappanames[EVALDISTR_Q], C->kappanames[EVALDISTR_DIM]);
	cov->q[size_idx] = nrow[EVALDISTR_Q] * ncol[EVALDISTR_Q] / dim;
	zaehler++;
      } 
      if (!PisNULL(EVALDISTR_P)) { // q
	if (ncol[EVALDISTR_P] != 1) 
	  SERR1("'%.50s' must be a vector", C->kappanames[EVALDISTR_P]);
	cov->q[size_idx] = nrow[EVALDISTR_P]  * dim;
 	zaehler++;
      } 
    }// kein else wegen zaehler !!

    if (!PisNULL(EVALDISTR_N)) { // r      
      cov->q[size_idx] = P0(EVALDISTR_N) * dim;
      zaehler++;
    } 
    if (zaehler != 1) SERR("exactly one of the parameters must be given");

  }
 
  //  if (!isRandom(sub)) SERR1("'%.50s' is not a distribution", NICK(sub));

  // PMI(sub);
  
  if ((err = CHECK_R(sub, dim)) != NOERROR) RETURN_ERR(err);
  assert(isnowRandom(sub));
 

  setbackward(cov, sub);  

  RETURN_NOERROR;
}

// S TRUCT(!isCovariance(cov) ? NULL : &(cov->key));



void range_EvalDistr(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[EVALDISTR_X] = RF_NEGINF;
  range->max[EVALDISTR_X] = RF_INF;
  range->pmin[EVALDISTR_X] = - 1e8;
  range->pmax[EVALDISTR_X] = 1e8;
  range->openmin[EVALDISTR_X] = true;
  range->openmax[EVALDISTR_X] = true;

  range->min[EVALDISTR_Q] = RF_NEGINF;
  range->max[EVALDISTR_Q] = RF_INF;
  range->pmin[EVALDISTR_Q] = - 1e8;
  range->pmax[EVALDISTR_Q] = 1e8;
  range->openmin[EVALDISTR_Q] = true;
  range->openmax[EVALDISTR_Q] = true;

  range->min[EVALDISTR_P] = 0;
  range->max[EVALDISTR_P] = 1;
  range->pmin[EVALDISTR_P] = 0;
  range->pmax[EVALDISTR_P] = 1;
  range->openmin[EVALDISTR_P] = false;
  range->openmax[EVALDISTR_P] = false;

  range->min[EVALDISTR_N] = 1;
  range->max[EVALDISTR_N] = RF_INF;
  range->pmin[EVALDISTR_N] = 1;
  range->pmax[EVALDISTR_N] = 1e8;
  range->openmin[EVALDISTR_N] = false;
  range->openmax[EVALDISTR_N] = false;

  range->min[EVALDISTR_DIM] = 1;
  range->max[EVALDISTR_DIM] = RF_INF;
  range->pmin[EVALDISTR_DIM] = 1;
  range->pmax[EVALDISTR_DIM] = 10;
  range->openmin[EVALDISTR_DIM] = false;
  range->openmax[EVALDISTR_DIM] = true;
}


int struct_EvalDistr(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  model *next = cov->sub[0],
    *sub = next;
  int  err,
    dim = ANYDIM;
   
  //  cov->simu.active = next->simu.active = false; 
  if (PL >= PL_DETAILS) { PRINTF("Struct EvalDistr\n"); }

  if ((err = STRUCT(sub, NULL)) != NOERROR) { RETURN_ERR(err); }
  if (PL >= PL_DETAILS) { PRINTF("Checking EvalDistr\n"); }
  
  if ((err = CHECK_R(sub, dim)) != NOERROR) RETURN_ERR(err);
    
  if (PL >= PL_DETAILS) {
    PRINTF("\n\nStruct EvalDistr (%s, #=%d), after 2nd check:",
	   NICK(sub), SUBNR); 
  }
  
  NEW_STORAGE(gen);

  if ((err = INIT(sub, 0, cov->Sgen)) != NOERROR) {
    RETURN_ERR(err); 
  }

  if (cov->rf == NULL) {
    int size = cov->q[0];
    if (cov->qlen > 1) size *= cov->q[1];

    if ((cov->rf = (double*) MALLOC(sizeof(double) * size)) 
	== NULL) RETURN_ERR(ERRORMEMORYALLOCATION); 
    cov->fieldreturn = wahr;
    cov->origrf = true;
  }
  
  //  cov->simu.active = sub->simu.active = true; 
  RETURN_NOERROR;
}


//void do_EvalDistr(model *cov, gen_storage *s){
//  assert(false);
//}


void Dummy(double VARIABLE_IS_NOT_USED *x,  INFO,
	   model VARIABLE_IS_NOT_USED *cov,
	   double VARIABLE_IS_NOT_USED *v) {
  BUG; //ERR("unexpected call of dummy function");
}


int check_dummy(model *cov) {
  model *next = cov->sub[0];
  defn *C = DefList + NEXTNR;
  int
    dim = Loctsdim(cov),
    err = NOERROR;
  Types
    nexttype = SYSTYPE(C->systems[0], 0);

  if (nexttype == RandomType) {
    if ((err = CHECK_R(next, dim) != NOERROR) && dim > 1)
      err = CHECK_R(next, 1);
    if (err != NOERROR) RETURN_ERR(err);
  } else {

#define nTF 3
#define ST ShapeType
#define ET EvaluationType
    Types
      types[nTF] = {ST, ST,     ProcessType},
      frames[nTF] ={ET, ET,     InterfaceType};
    domain_type 
      dom[nTF] ={XONLY, KERNEL, XONLY};
    int 
      start = 1,
      end = nTF;
    domain_type
      prevdom = PREVDOM(0),
      nextdom = DOM(C->systems[0], 0);
    isotropy_type iso  = OWNISO(0);
      //CoordinateSystemOf(PREVISO(0));
    if (equalsUnreduced(iso)) iso = CoordinateSystemOf(iso);
    
    switch(nexttype) {  
    case TcfType: case PosDefType: case VariogramType: case NegDefType:
    case PointShapeType : case ShapeType:  case TrendType:
    case RandomOrShapeType :case MathDefType:
      end = 2; break;
    case ProcessType:  case NormedProcessType:  case BrMethodType:
    case SmithType:  case SchlatherType:  case PoissonType:
    case PoissonGaussType:  case GaussMethodType:
      start = 3; break;
    default :
      if (nexttype != ManifoldType && nexttype != OtherType)
	// e.g. InterfaceType
	BUG;
    }

    if (equalsKernel(nextdom) ||
	(isNegDef(nexttype) && isAnySpherical(iso))) {
      if (start==1) start = 2;
      else BUG;
    } else if (equalsXonly(nextdom)) {
      if (end == 2) end = 1;
    }
   
    ASSERT_LOC_GIVEN;
    for (int t=start-1; t<end; t++) {
      // printf("f=%d (%s, %s, %s)\n", t,	     TYPE_NAMES[frames[t]], TYPE_NAMES[types[t]],	     DOMAIN_NAMES[dom[t]]);
      if ((err = CHECK(next, dim, dim, types[t], dom[t], iso, SUBMODEL_DEP,
		       frames[t]))
	  == NOERROR) break;
      //      PMI1(cov);
      // printf("end = %d\n", end);
    }
    if (err != NOERROR) RETURN_ERR(err);
  }
  
  setbackward(cov, next);  
  VDIM0 = next->vdim[0]; 
  VDIM1 = next->vdim[1]; 

  RETURN_NOERROR;
}


int struct_dummy(model VARIABLE_IS_NOT_USED *cov, model VARIABLE_IS_NOT_USED **newmodel){    
  RETURN_NOERROR;
}


#define RFGET_UP 0
#define RFGET_REGISTER 1
#define RFGET_SUB 0
void RFget(double VARIABLE_IS_NOT_USED *x, INFO, model *cov, double *v){
  getStorage(s ,   get);
  GETSTOMODEL;
  model *get_cov = STOMODEL->get_cov;
  int i,
    nr = MODELNR(STOMODEL->get_cov),
    param_nr = s->param_nr,
    *idx = s->idx,
    size = s->size;

  if (DefList[nr].kappatype[param_nr] == REALSXP) {
    double *p = PARAM(get_cov, param_nr);    
    if (s->all) {
      for (i=0; i<size; i++) v[i] = p[i];
    } else {
      for (i=0; i<size; i++) v[i] = p[(int) idx[i]];
    }
  } else if (DefList[nr].kappatype[param_nr] == INTSXP) {    
    int *p = PARAMINT(get_cov, param_nr);
    if (s->all) {
      for (i=0; i<size; i++) v[i] = (double) p[i];
    } else {
      for (i=0; i<size; i++) v[i] = (double) p[idx[i]];
    }
  } else BUG;
}



int SearchParam(model *cov, get_storage *s, model_storage *STOMODEL) {
  model *orig;
  int i, subcov,
    up = P0INT(RFGET_UP);
  if (PisNULL(RFGET_REGISTER)) {
    orig = cov->calling;
    if (orig == NULL) SERR("register not given");
    for (i=0; i<up; i++) {
      orig = orig->calling;
      while (orig != NULL && orig->user_given == ug_internal)
	orig = orig->calling;
      if (orig == NULL) SERR1("value of '%.50s' too high", KNAME(RFGET_UP));
    }
  } else {
    int nr = P0INT(RFGET_REGISTER);
    if (nr<0 || nr>MODEL_MAX) SERR("invalid register number");
    if (up != 0) SERR1("'%.50s' may not be given.", KNAME(RFGET_UP));
    orig = KEY()[nr];
  }
  STOMODEL->orig = orig;
  
  while (true) {
    while (DefList[MODELNR(orig)].maxsub > 0 && orig != NULL &&
	   orig->user_given == ug_internal) 
      orig =  COVMODELKEYS_GIVEN(orig) &&
	(MODELNR(orig) == PLUS || MODELNR(orig)==MULT || MODELNR(orig)==MPPPLUS)
	? orig->Smodel->keys[0]
	: orig->key != NULL ? orig->key
	: orig->sub[0];
    if (orig == NULL || DefList[MODELNR(orig)].maxsub == 0) 
      SERR("model does not match the request or is too complicated");
    if (COVNR != MODELNR(orig)) 
      SERR2("(sub)models (%.50s, %.50s) do not match",
	    NICK(cov), NICK(orig));
    int nsub = DefList[MODELNR(orig)].maxsub;
    for (subcov=0; subcov < nsub; subcov++) 
      if (cov->sub[subcov] != NULL) break;   
    if (subcov < nsub) { // submodel found
      if (orig->sub[subcov] == NULL) 
	SERR2("(sub)models (%.50s; %.50s) do not match",
	      NICK(cov), NICK(orig));
      cov  = cov->sub[subcov];
      orig = orig->sub[subcov];
    } else {
      int
	kappas = DefList[MODELNR(orig)].kappas;
      for (i=0; i < kappas; i++) 
	if (cov->kappasub[i] != NULL) break;         
      if (i < kappas) { // param submodel found
	if (orig->kappasub[i] == NULL) 
	  SERR2("parameter models of %.50s and %.50s do not match",
		NICK(cov), NICK(orig));
        cov  = cov->kappasub[i];
	orig = orig->kappasub[i];   
      } else {
	for (i=0; i < kappas; i++) 
	  if (cov->kappasub[i] != NULL) break;         
	if (i >= kappas) SERR("no parameter given");
	defn *C = DefList + COVNR;
	if (C->kappatype[i] == REALSXP) s->all = P(i)[0] == 0;
	else if (C->kappatype[i] == INTSXP) s->all = PINT(i)[0] == 0;
	else SERR("only numeric parameters are allowed");
	if (s->all) {
	  s->vdim[0] = orig->nrow[i];
	  s->vdim[1] = orig->ncol[i];
	  s->size = s->vdim[0] * s->vdim[1];
	} else {
	  int k;
	  s->size = s->vdim[0] = cov->nrow[i];
	  s->vdim[1] = cov->ncol[i];
	  if (cov->ncol[i] != 1) SERR("only vectors of indices are allowed");
	  assert(s->idx == NULL);
	  s->idx = (int *) MALLOC(sizeof(int) * s->size);
	  if (C->kappatype[i] == REALSXP)
	    for (k=0; k<s->size; k++) s->idx[k] = ((int) P(i)[k]) - 1;
	  else 
	    for (k=0; k<s->size; k++) s->idx[k] = PINT(i)[k] - 1;
	}
	STOMODEL->get_cov = orig;
	s->param_nr = i;
	RETURN_NOERROR;
      }
    }
  } // while true
  BUG;
  RETURN_ERR(ERRORFAILED); // nur fuer compiler
}
  
int check_RFget(model *cov) {

  BUG; /// todo:  Code ist noch nicht ausgegoren !!

  //model *orig, *sub;
  int i, err;
  //    len = ((idx != NULL) ? cov->nrow[RFGET_IDX]
  //	  : get->get_cov->ncol[param_nr] * get->get_cov->nrow[param_nr]);
  if (cov->Sget != NULL) RETURN_NOERROR;

  kdefault(cov, RFGET_UP, 0);
  NEW_STORAGE(get); // fuehrt zu fehler i.A.
  NEWSTOMODEL;
  getStorage(s ,   get);
  GETSTOMODEL;

  if ((err = SearchParam(cov, s, STOMODEL)) != NOERROR) RETURN_ERR(err);
  for (i=0; i<2; i++) cov->vdim[i] = s->vdim[i];
  RETURN_NOERROR;
}
  

void range_RFget(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[RFGET_UP] = 0;
  range->max[RFGET_UP] = RF_INF;
  range->pmin[RFGET_UP] = 0;
  range->pmax[RFGET_UP] = 100;
  range->openmin[RFGET_UP] = false;
  range->openmax[RFGET_UP] = true;

  range->min[RFGET_REGISTER] = 0;
  range->max[RFGET_REGISTER] = MODEL_MAX;
  range->pmin[RFGET_REGISTER] = 0;
  range->pmax[RFGET_REGISTER] = MODEL_USER;
  range->openmin[RFGET_REGISTER] = false;
  range->openmax[RFGET_REGISTER] = false;
}


int struct_RFget(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  int i, err;
  
  NEW_STORAGE(get); // Fehler ?!
  NEWSTOMODEL;
  getStorage(s ,   get); 
  GETSTOMODEL;
  if ((err = SearchParam(cov, s, STOMODEL)) != NOERROR) RETURN_ERR(err);
  for (i=0; i<2; i++) {
    if (cov->vdim[i] != s->vdim[i])
      SERR("mismatch of dimensions when constructing the model");
  }
  cov->fieldreturn = wahr; // seltsam !
  cov->origrf = false;  
  RETURN_NOERROR;
}




//void do_Rfget(model *cov, gen_storage *s){
//  assert(false);
//}






//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

model *get_around_gauss(model *cov) {
  model *next= cov;
  if (NEXTNR == SCHLATHERPROC) next = next->sub[0]; // kein else 
  if (NEXTNR == GAUSSPROC) next = next->sub[0];

  if (isGaussMethod(next) || equalsBernoulliProcess(next)) {
    if (NEXTNR == AVERAGE_USER){
      if (next->sub[0] == NULL) ERR("covariance cannot be calculated (yet) for arbitrary shape functions.");
      next = next->sub[next->sub[0] == NULL];
      if (NEXTNR == AVERAGE_INTERN) next = next->sub[next->sub[0] == NULL];
    } 
    else if (NEXTNR == CE_CUTOFFPROC_USER) {
      next = next->sub[0];      
      if (NEXTNR == CE_CUTOFFPROC_INTERN) next = next->sub[0];
    }
    else if (NEXTNR == CE_INTRINPROC_USER) {
      next = next->sub[0];   
      if (NEXTNR == CE_INTRINPROC_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == NUGGET_USER) {
      next = next->sub[0];   
      if (NEXTNR == NUGGET_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == HYPERPLANE_USER) {
      next = next->sub[0];   
      if (NEXTNR == HYPERPLANE_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == SPECTRAL_PROC_USER) {
      next = next->sub[0];   
      if (NEXTNR == SPECTRAL_PROC_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == TBM_PROC_USER) {
      next = next->sub[0];   
      if (NEXTNR == TBM_PROC_INTERN) next = next->sub[0];   
    }
    else if (NEXTNR == RANDOMCOIN_USER) {
      if (next->sub[0] == NULL) ERR("covariance cannot be calculated (yet) for arbitrary shape functions.");
      next = next->sub[next->sub[0] == NULL];   
      if (NEXTNR == AVERAGE_INTERN) next = next->sub[next->sub[0] == NULL];
    } else {
      BUG;
    }
  }
  return next;
}



model *get_around_max_stable(model *cov) {
  model *next = cov;

  if (isBrMethod(next)) {
    next = next->sub[0];
    if (MODELNR(next->calling) == BROWNRESNICKPROC && isBrMethod(next)) {
      next = next->sub[0];
    }
  } 
  return next;
}



int check_fctn_intern(model *cov, Types type, bool close,
		     usekernel_type kernel, 
		     int rows, int cols, Types frame) {
  //coordinates define the return values, i.e. always LocFoo instead of LocFooY 
  model
    *next = cov->sub[0],
    *sub = cov->key == NULL ? next : cov->key;
  ASSERT_LOC_GIVEN;

  int err = ERRORFAILED,   // +  (int) nr_coord_sys_proj], 
    dim =  Loctsdim(cov);
  assert(dim == OWNLOGDIM(0));
  
  domain_type
    firstdomain = kernel==forceKernel ? KERNEL : XONLY,
    lastdomain = kernel==forceKernel ? KERNEL : XONLY; 
  if (kernel==allowKernel) {
    if (isNegDef(type) && isAnySpherical(OWNISO(0)))firstdomain = KERNEL;
    if (!(isTrend(type) || isProcess(type))) lastdomain = KERNEL;
  }

  // PMI(cov)
  // firstdomain = XONLY;
  // if (LocHasY(cov)) firstdomain = lastdomain = KERNEL;
  
  //   printf("firstdomain %d %d %d\n", firstdomain, lastdomain, kernel);

  
  // assert(firstdomain == XONLY);

  /*
    int endfor = 0,
    isotropy_type iso[4];
    iso[endfor++] = type == ShapeType 
    ? CoordinateSystemOf(OWNISO(0)) 
    : S ymmetricOf(OWNISO(0));
  
    if (isIsoMismatch(iso[endfor-1])) BUG;
    int end_frame = sub = = next;
    for (int i=0; i < endfor; i++) {
    for (k=XO NLY; k<=lastdomain; k++) {
    Types frame = V ariogramType;
    for (int f=0; f<=end_frame; f++) {
    if ((err = CHECK(sub, dim, OWNXDIM(0), type, 
    (domain_type) k, iso[i], SUBMODEL_DEP, frame))
    //sub!=next || isVariogram(sub) ? V ariogramType : Any Type))
    == NOERROR) break;
    if (err == NOERROR) break;
    frame = Any Type;
    }
    if (err == NOERROR) break;
    }
    if (err == NOERROR) break;
    }
  */

 
  //  Types frame = isProcess(type) ? EvaluationType : EvaluationType;
  //  int end_frame = sub == next && equalsEvaluation(frame);

  //  assert(frame = EvaluationType);
  //  printf("hrt\n");

  isotropy_type iso = SymmetricOf(OWNISO(0));
  for (int j=0; j<=1; j++) {
    if (equalsIsoMismatch(iso)) BUG;
    for (int k=firstdomain; k<=lastdomain; k++) {
      //printf("j=%d k=%d %d %d\n",j, k, firstdomain, lastdomain);
      // for (int f=0; f<=end_frame; f++) {
      if ((err = CHECK(sub, dim, OWNXDIM(0),
		       type, 
		       (domain_type) k, iso, SUBMODEL_DEP, frame))
	  //sub!=next || isVariogram(sub) ? V ariogramType : Any Type))
	  == NOERROR) break;
      //PMI1(cov);
    }
    iso = CoordinateSystemOf(OWNISO(0));
  }

  if (err != NOERROR) RETURN_ERR(err);
  setbackward(cov, sub); 

  // memory set according to the submodel as this model will
  // be evaluated (cov clear; fctn as function; predict: cov-matrix!)
  if ((err = alloc_fctn(cov, dim, VDIM0 * VDIM1)) != NOERROR)
    RETURN_ERR(err);

  // this is how cov will forward the result
  if (rows > 0) VDIM0 = rows;
  if (cols > 0) VDIM1 = cols;

  if (sub->pref[Nothing] == PREF_NONE) SERR("given model cannot be evaluated");
  
  if (cov->q == NULL) {
    int d,
      len=1; // # of "simulations" (repetitions!)
    bool grid = Locgrid(cov);
    coord_type xgr = Locxgr(cov);
    if (grid) len += dim; else len ++;      
    for (int i=0; i<2; i++) len += (int) (cov->vdim[i] > 1);
    QALLOC(len);
    
#define VDIMS								\
    for (int i=0; i<2; i++) if (cov->vdim[i] > 1) cov->q[d++] = cov->vdim[i]
#define LOCS if (grid) {						\
      for (int i=0; i<dim; i++) cov->q[d++] = xgr[i][XLENGTH];		\
    } else {								\
      cov->q[d++] = Loctotalpoints(cov);				\
    }      
    
    d = 0;
    if (close) {
      VDIMS;
      LOCS;	
    } else {
      LOCS;
      VDIMS;
    }
    cov->q[d] = 1;
    assert(d == len-1);
  }
  
  RETURN_NOERROR;
}


void Cov(double VARIABLE_IS_NOT_USED *x, INFO, model *cov, double *v){
  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  Covariance(sub, v);
  //  for (int i=0; i<512; i++) printf("%f ", v[i]); assert(false);
  //APMI(cov);
}

isotropy_type which_system[nr_coord_sys] =
  { CARTESIAN_COORD, ISO_MISMATCH, CARTESIAN_COORD, EARTH_COORD,
    SPHERICAL_COORD, CARTESIAN_COORD, CARTESIAN_COORD, ISO_MISMATCH };

int check_cov_intern(model *cov, Types type, bool close,
		     usekernel_type kernel) {
  globalparam *global = &(cov->base->global);
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;


  // printf("ispRocess %d\n", isProcess(sub));
  
  if (isProcess(sub)) {
    int err,
      dim =  Loctsdim(cov);
    
    err = CHECK_THROUGHOUT(sub, cov, ProcessType, XONLY,
			   which_system[global->coords.coord_system],
			   SUBMODEL_DEP,
			   EvaluationType);
    //  int err = CHECK(sub, cov->tsdim, cov->xdimown, ProcessType, XONLY, 
    //		    OWNISO(0), SUBMODEL_DEP, 
    //		    cov->frame = = Any Type ? Any Type : GaussMethodType);
    if (err != NOERROR) RETURN_ERR(err);
    setbackward(cov, sub);
    VDIM0 = sub->vdim[0];
    VDIM1 = sub->vdim[1];
    if ((err = alloc_fctn(cov, dim, VDIM0 * VDIM1)) != NOERROR) RETURN_ERR(err);
    RETURN_NOERROR;
  } else return check_fctn_intern(cov, type, close, kernel, 0,0,EvaluationType);
}

int check_cov(model *cov) {
  globalparam *global = &(cov->base->global);
  return check_cov_intern(cov, PosDefType, global->general.vdim_close_together,
			  LocHasY(cov) ? forceKernel : noKernel);
}

int struct_cov(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  assert(cov->sub[0] != NULL);
  model
    *next = cov->sub[0],
    *sub = get_around_gauss(next);

  // fehlt : Poisson

  if (sub != next) {
    if ((err = COVNR == COVMATRIX ? check_covmatrix(cov) : check_cov(cov))
	!= NOERROR) RETURN_ERR(err);

    //    if (cov->Seval == NULL) { PMI0(cov); crash(); }
    
    assert(cov->Sfctn != NULL);
    
    ONCE_NEW_STORAGE(gen);

    err = INIT(next, 0, cov->Sgen);
    RETURN_ERR(err);
  }
  RETURN_NOERROR;
}


int init_cov(model *cov, gen_storage *s) {
  // darf nur von Likelihood aus aufgerufen werden
  if (hasAnyEvaluationFrame(cov)) {
    BUG;
    assert(cov->key == NULL);
    return INIT(cov->sub[0], 0, s);
  }

  RETURN_ERR(ERRORFAILED);
}

void CovMatrix(double VARIABLE_IS_NOT_USED *x, INFO, model *cov, double *v){
  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  if (DefList[SUBNR].is_covmatrix(sub)) DefList[SUBNR].covmatrix(sub, true, v);
  else StandardCovMatrix(sub, true, v);
}

int check_covmatrix(model *cov) {
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key;
  assertNoLocY(cov);
  ASSERT_LOC_GIVEN;
  int err = NOERROR,
    rows, cols,
    dim = Loctsdim(cov), // !! not GettLoctsdim
    total = Loctotalpoints(cov); // !! not Loctotalpoints
  isotropy_type iso = SymmetricOf(PREVISO(0));

  assert(dim == OWNLOGDIM(0));
  if (LocDist(cov) && LocxdimOZ(cov) == 1) {
    iso = PREVISO(0);
    iso = isCartesian(iso) ? ISOTROPIC
      : isEarth(iso) ? EARTH_ISOTROPIC
      : isSpherical(iso) ? SPHERICAL_SYMMETRIC
      : ISO_MISMATCH;
  } else {
    if (PREVXDIM(0) != PREVLOGDIM(0)) BUG;
  }

  if ((err = CHECK(sub, dim, OWNXDIM(0), PosDefType, XONLY, iso,
		   SUBMODEL_DEP, EvaluationType)) != NOERROR) {
    if ((err = CHECK(sub, dim, OWNXDIM(0), 
		     PosDefType, KERNEL, CARTESIAN_COORD, SUBMODEL_DEP,
		     EvaluationType)) != NOERROR) {
      //      APMI(cov)
      if ((err = CHECK(sub, dim, OWNXDIM(0), VariogramType, XONLY,
		       // iso,
		       SymmetricOf(PREVISO(0)),
		       SUBMODEL_DEP, EvaluationType)) != NOERROR) {
	RETURN_ERR(err);
      }
    }
  }

  setbackward(cov, sub);  
  rows = VDIM0 = sub->vdim[0]; 
  cols = VDIM1 = sub->vdim[1]; 

  if (cov->q == NULL) {
    QALLOC(2);
    cov->q[0] = total * rows;
    cov->q[1] = total * cols;
  }

  if ((err = alloc_fctn(cov, dim, rows * cols)) != NOERROR) RETURN_ERR(err);
  
  RETURN_NOERROR;
}


void Pseudomadogram(double VARIABLE_IS_NOT_USED *x, INFO, model *cov,
		    double *v){
  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  Pseudomadogram(cov->key == NULL ? cov->sub[0] : cov->key,
		 P0(PSEUDO_ALPHA), v);
}


int check_pseudomado(model *cov) {
  globalparam *global = &(cov->base->global);
  kdefault(cov, PSEUDO_ALPHA, 2.0);
  double alpha = P0(PSEUDO_ALPHA);
  return check_cov_intern(cov,
			  alpha == VARIOGRAM ? VariogramType :
			  alpha == COVARIANCE ? PosDefType :
			  VariogramType, // NegDefType ??? 22.1.21
			  global->general.vdim_close_together,
			  LocHasY(cov) ? forceKernel : noKernel);
}


void range_pseudomado(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[PSEUDO_ALPHA] = 2.0;
  range->max[PSEUDO_ALPHA] = 2.0;
  range->pmin[PSEUDO_ALPHA] = 2.0;
  range->pmax[PSEUDO_ALPHA] = 2.0;
  range->openmin[PSEUDO_ALPHA] = false;
  range->openmax[PSEUDO_ALPHA] = false;
}


int check_pseudovario(model *cov) {
  globalparam *global = &(cov->base->global);
  kdefault(cov, PSEUDO_ALPHA, 2.0);
  return check_cov_intern(cov, VariogramType,
			  global->general.vdim_close_together,
			 LocHasY(cov) ? forceKernel : noKernel );
}

void range_pseudovario(model VARIABLE_IS_NOT_USED *cov, range_type* range){
  range->min[PSEUDO_ALPHA] = 0.0;
  range->max[PSEUDO_ALPHA] = 2.0;
  range->pmin[PSEUDO_ALPHA] = 0.1;
  range->pmax[PSEUDO_ALPHA] = 2.0;
  range->openmin[PSEUDO_ALPHA] = true;
  range->openmax[PSEUDO_ALPHA] = false;
}


void Variogram(double VARIABLE_IS_NOT_USED *x, INFO, model *cov, double *v){
  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array  
  Variogram(cov->key == NULL ? cov->sub[0] : cov->key, v);
}


int check_vario(model *cov) {
  globalparam *global = &(cov->base->global);
//  printf("lochasy = %d\n",  LocHasY(cov) );
  return check_cov_intern(cov, VariogramType,
			  global->general.vdim_close_together,
			  LocHasY(cov) ? forceKernel : noKernel);
//  if (err == NOERROR) PMI(cov) else M ERR(err);
//    return er r;
}


int struct_variogram(model *cov, model VARIABLE_IS_NOT_USED **newmodel){
  int err;
  model *sub,
    *next = cov->sub[0];
 
  if ((sub = get_around_max_stable(next)) == next) sub = get_around_gauss(next);
  if (sub != next) {
    if ((err = covcpy(&(cov->key), sub)) != NOERROR) RETURN_ERR(err);       
    sub = cov->key;
    SET_CALLING(sub, cov);

    if ((err = CHECK(sub, Loctsdim(cov), OWNXDIM(0), VariogramType,
		     NEXTDOM(0), NEXTISO(0), cov->vdim,
		     EvaluationType)) != NOERROR) {
      RETURN_ERR(err);
    }
  }
  
  if (!isnowVariogram(sub))
    SERR(sub != next ? "variogram model cannot be determined"
	 :"not a variogram model");
  
  RETURN_NOERROR;
} 

   
// bool close = gl obal->general.vdim_close_together;


void Fctn(model *cov, bool ignore_y, double *v) {
  model *sub = cov->sub[0];
  if (sub == NULL) { BUG;}
  FctnIntern(cov, cov, sub, ignore_y, v);
}


void Fctn(double VARIABLE_IS_NOT_USED *x, INFO, model *cov, double *v) {
  Fctn(cov, false, v); // user !
}


void FctnIntern(model *cov, model *covVdim, model *genuine, bool ignore_y,
		double *v){
  assert(cov != NULL);
  globalparam *global = &(cov->base->global);
  if (v==NULL) return; // EvaluateModel needs information about size
  //                      of result array
  //  PMI0(cov); // PMI0(genuine->calling);  PMI0(genuine);
   assert(hasLikelihoodFrame(cov) || hasInterfaceFrame(cov) || isProcess(cov));

   // assert(XDIM(PREVSYSOF(genuine), 0) == LOGDIM(PREVSYSOF(genuine), 0));
   if (XDIM(PREVSYSOF(genuine), 0) != LOGDIM(PREVSYSOF(genuine), 0)) BUG;

  
  assert(genuine != NULL);
  //  PMI(cov);
  // crash();
  
  // assert(cov->key == NULL);
  bool kernel = isKernel(PREVSYSOF(genuine));
  FINISH_START(cov, covVdim, !kernel && !ignore_y, ignore_y, 1, 0, 1, 0);
  info[INFO_EXTRA_DATA_X] = !kernel && !ignore_y;
  info[INFO_EXTRA_DATA_Y] = !ignore_y;
 
  
  assert(fctn != NULL);
  double 
    *cross = fctn->cross,
    *zero = ZERO(covVdim);

  y = ygiven ? fctn->y : ZERO(cov);
  Types frame = genuine->frame;
  if (hasAnyEvaluationFrame(genuine)) genuine->frame = ShapeType;

 if (vdimSq > fctn->rowscols) {PMI(cov); BUG;} //

  assert(hasFullXdim(ISO(PREVSYSOF(genuine), 0))
	 ||
	 (ISO(PREVSYSOF(genuine), 0) == ISO(SYSOF(genuine), 0)
	  &&  isTrend(SYSTYPE(PREVSYSOF(genuine), 0) )
	  ) 
	 );


  // printf("\n%d \n", totX);
  
#define FUNCTION FCTN(x, info, genuine, v)
#define FUNCTION2 FCTN(x, info, genuine, cross); VDIM_LOOP(cross[u])
#define FUNCTION_Y NONSTATCOV(x, y, info, genuine, v);
#define FUNCTION2_Y NONSTATCOV(x, y, info, genuine, cross); VDIM_LOOP(cross[u])
   
  PERFORM_PREPARE;

  //ygiven = 1 kernel=0 grid=0 0 kernel.genuine=0
  // printf("ygiven = %d %d kernel=%d grid=%d %d kernel.genuine=%d\n", ygiven, ignore_y,  kernel, grid, gridY, isKernel(PREVSYSOF(genuine)));
  // PMI0(genuine);
  
  PERFORM(FUNCTION, FUNCTION2, FUNCTION_Y, FUNCTION2_Y);
  STANDARD_ENDE;

  genuine->frame = frame; // ja nicht loeschen!!
}


void FctnExtern(model *cov, model *genuine, bool ignore_y, double *v){
  Types frame = cov->frame;
  int dim =  Loctsdim(cov);
  if (alloc_fctn(cov, dim, cov->vdim[0] * cov->vdim[1]) != NOERROR)
    XERR(ERRORMEMORYALLOCATION);
  cov->frame = LikelihoodType; // dummy, just to please next function
  FctnIntern(cov, cov, genuine, ignore_y, v);
  cov->frame = frame;
  fctn_DELETE(&(cov->Sfctn));
}

 
int check_fctn_intern(model *cov) {
  assert(COVNR == SHAPE_FCT_PROC || COVNR == PROD_PROC || COVNR == FCTNFCTN);
  globalparam *global = &(cov->base->global);
#define nTypesTF 2  
  Types T[nTypesTF] = {TrendType, ShapeType},
    F[nTypesTF] = {TrendType, LikelihoodType};
  int i, err;
  usekernel_type kernel = COVNR == FCTNFCTN && LocHasY(cov)
    ? forceKernel : noKernel;

  //printf("checkfctn %d %d %d\n", COVNR == FCTNFCTN, LocHasY(cov), kernel);
  
  for (i=0; i<nTypesTF; i++) {
    err = check_fctn_intern(cov, T[i], global->general.vdim_close_together, 
			   kernel, 0, 0, F[i]);
    if (err == NOERROR) RETURN_ERR(err);
  }
  // BUG vor 3.4.2018 -- warum?
  ONCE_EXTRA_STORAGE;
  RETURN_ERR(err);
}

 
int check_fctn(model *cov) {
  cov->base->rawConcerns = onlyOneDataSetAllowed;
  return check_fctn_intern(cov);
}



/* ****************************** */
/*             PREDICT            */
/* ****************************** */


void kappapredict(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc) {
  if (i <= LIKELIHOOD_LAST) {
    kappalikelihood(i, cov, nr, nc);
    return;
  }
  // if (i==PREDICT_GIVENIDX || i==PREDICT_PREDICTIDX) {
  //    *nr = SIZE_NOT_DETERMINED;
  //    *nc = 1;
  //  } else
  *nr = *nc = (i==PREDICT_GIVEN) ? SIZE_NOT_DETERMINED : -1;
}

void predict(double *N, int *idx, model *cov, double *v) {
  // printf("predict\n");
  /*
    N > 0 : dim(N) == 2 && N[0] = length(idx for x);
                           N[1] = length(idx for data)  
           *idx = c(idx for x, idx for data)
    N == 0.0 : idx[0] = set number
   */

  // PMI1(cov);
  
  model *sub = cov->key == NULL? cov->sub[PREDICT_CONDITIONING] : cov->key ;
 

  if (v==NULL) {
    GETSTORAGE(L ,  sub,   likelihood); 
    assert(L != NULL);
    cov->base->set = 0;
    listoftype *datasets = L->datasets;
    int
      vdim = VDIM0,
      ncol =  NCOL_OUT_OF(datasets),
      repet = ncol / vdim;
    assert(cov->qlen == 3 && cov->q != NULL);
    cov->q[cov->qlen - 1] = repet;
    return; // EvaluateModel needs information about size
   //                      of result array, given in cov->q
  }
 
  int n = (int) *N;

  if (n == 0) {
    if (idx[0] > LocSets(cov)) BUG;
    // printf("idx=%d\n", idx[0]);
    cov->base->set = idx[0] - 1;// idx[0] has numbering of R
    assert(cov->base->set >= 0);
    if (SUBNR == GAUSSPROC) {
      gauss_predict(cov, v);
      cov->base->set = 0;
      return;
    }
  }

  //  else if (n > 0) { // not used yet
  //    PFREE(PREDICT_IDX);
  //    PALLOC(PREDICT_IDX, n, 1);
  //    BUG;
  //    // idx hat die R indizierung 
  //    // Auf P(PREDICT_IDX) kopieren und um 1 reduzieren
  //    return;
  //  } 

  
  BUG;
}

int check_predict(model *cov) {
  assert(PREDICT_CONDITIONING == 0);
  if (!LocHasY(cov)) loc_set(PVEC(PREDICT_GIVEN) , LocP(cov));
  if (cov->q == NULL) {
    QALLOC(3); // wichtig for dem checke_likeli
    cov->q[0] = Loctotalpoints(cov);
  }
  int err = check_likelihood(cov);
  if (err != NOERROR) RETURN_ERR(err);
 
  assert(cov->prevloc != NULL);
  bool distances = LocDist(cov); // ehem Prev
  model *sub = cov->sub[PREDICT_PREDICT];
  int dim = Loctsdim(cov);
  domain_type dom = distances && LocxdimOZ(cov) == 1 ? XONLY : KERNEL;
  Types frame = LikelihoodType;
    // SUBNR == PLUS || isProcess(sub) ? LikelihoodType : EvaluationType;
  
  if (isProcess(sub)) {
    err = CHECK(sub, dim, dim, ProcessType, dom,
		distances ? ISOTROPIC : UNREDUCED,
		cov->vdim, frame);
  } else {
    err = CHECK(sub, dim, dim, PosDefType, dom,
		distances ? ISOTROPIC : CoordinateSystemOf(PREVISO(0)),
		cov->vdim,
		frame);
    if (err != NOERROR) {
      err = CHECK(sub, dim, dim, NegDefType, dom,
		  distances ? ISOTROPIC : SymmetricOf(PREVISO(0)),
		  cov->vdim,
		  frame);
    }
  }
  if (err != NOERROR) RETURN_ERR(err);
  if (VDIM0 !=sub->vdim[0] || VDIM1 != sub->vdim[1])
    SERR("m-dimensionality of 'model' and 'err.model' do not match");
  
    
  if ((err = alloc_fctn(sub, dim, VDIM0 * VDIM1)) != NOERROR)
    RETURN_ERR(err);
  
  setbackward(cov, sub); 
  cov->q[1] = VDIM0;
 
  ONCE_EXTRA_STORAGE;
  RETURN_NOERROR;
}

  

int struct_predict(model *cov, model VARIABLE_IS_NOT_USED  **newmodel){
  int err = struct_likelihood(cov, newmodel);
  if (err != NOERROR) RETURN_ERR(err);
  // return struct_cov(cov, newmodel);
  RETURN_NOERROR;
}

void range_predict(model VARIABLE_IS_NOT_USED *predict, range_type* range){
  range_likelihood(predict, range); 
}
