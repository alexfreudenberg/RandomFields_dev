
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields -- init part

 Copyright (C) 2001 -- 2017 Martin Schlather

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



// to do: V3.1+: Randverteilungen der RPs in cov, D etc implementieren;
//        do -> random part; init wird innerhalb von do aufgerufen,
//        falls nicht initialisiert oder random submodels??
// to do: MLE: random parameters einsammeln


#include <stdio.h>  
#include <string.h>
#include "RF.h"
#include "operator.h"
#include "startGetNset.h"
#include "primitive.h"
#include "init.h"
#include "questions.h"
#include "shape.h"
#include "QMath.h"


int
  PLUS,
  USER,    
  FIRSTDOLLAR, LASTDOLLAR, 
  ISO2ISO0, SP2SP0, SP2ISO0, S2ISO0, S2SP0, S2S0, SId0, E2E0, E2EIso0, 
  E2Sph0, E2SphIso0, Sph2Sph0, Sph2SphIso0,  FIRSTGATTER0, LASTGATTER0,
  FIRST_PLANE, LAST_PLANE, EARTHKM2CART, EARTHMILES2CART,
  EARTHKM2GNOMONIC, EARTHMILES2GNOMONIC,
  EARTHKM2ORTHOGRAPHIC, EARTHMILES2ORTHOGRAPHIC,  
  FIRST_TRAFO, LAST_TRAFO, MATHDIV;

defn *DefList=NULL;
int currentNrCov=UNSET,
  gaussmethod[Forbidden+1];
double ONE = 1;
char CovNames[MAXNRCOVFCTS][MAXCHAR], CovNickNames[MAXNRCOVFCTS][MAXCHAR],
  *FREEVARIABLE= (char*) "...";


//                        CE         CO         CI        TBM       Sp
//                        di         sq       trend        av       n
//                        mpp         Hy       spf        any     forbidden
pref_type PREF_ALL = {PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, 
		      PREF_BEST, PREF_BEST, PREF_NONE, PREF_BEST, PREF_BEST,
		      PREF_BEST, PREF_BEST, PREF_NONE, // specific
		                                       PREF_BEST, PREF_BEST},
  PREF_NOTHING = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST, // nothing
		                                              PREF_NONE},
  PREF_TREND =  {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
		  PREF_NONE, PREF_NONE, PREF_BEST, PREF_NONE, PREF_NONE,
		  PREF_NONE, PREF_NONE, PREF_NONE, PREF_BEST, // nothing
		                                              PREF_NONE},
  PREF_AUX = {PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, 
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE,
	      PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE, PREF_NONE},
  PREF_MATHDEF = {PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST, PREF_BEST,
		  PREF_BEST, PREF_BEST, PREF_BEST, PREF_NONE, PREF_NONE,
		  PREF_BEST, PREF_BEST, PREF_NONE, PREF_BEST, PREF_NONE};
 

const char *CAT_TYPE_NAMES[OtherType + 1] = {
  // TcfType, PosDefType, VariogramType, NegDefType,
  // PointShapeType, ShapeType, TrendType, RandomOrShape, Manifold,
  // ProcessType, GaussMethodType, NormedProcessType, BrMethodType,
  // SmithType, SchlatherType, PoissonType, PoissonGaussType,
  // RandomType, InterfaceType, MathDefType, OtherType
  // not including BadType, ...
    "RM", "RM", "RM", "RM",
    "RM", "RM", "RM", "RM", "RM", 
    "RP", "RP", "RP", "RP",
    "RP", "RP", "RP", "RP",
    "RR",  "RF", "R.", "RO"},
  
  *REG_NAMES[MODEL_MAX+1] = {"reg0", "reg1", "reg2", "reg3", "reg4", 
			    "reg5", "reg6", "reg7", "reg8", "reg9",
			     "user", "cov", "covmatrix", "variogram",
			     "pseudovariogram",
			     "function", "distribution", "calculation",
			     "unused", "unused",
			     "auxiliary", "intern", "split", "gui", "mle",
			     "mlesplit", "mletrend", "mlebounds",
		 	    "kriging", "conditional", "error model"},
  *POSITIVITY_NAMES[(int) pt_mismatch + 1] = 
			      {"pt-wise positive definite", 
			       "pt-wise indefinite", 
			       "pt-wise negative definite", "pt-wise zero", 
			       "pt-wise param dependent", 
			       "pt-wise submodel dependent", 
			       "pt-wise undefined", "pt-wise unknown", 
			       "pt-wise mismatch"},
  **LIST_OF_NAMES[nNamesOfNames] = // see also NAME_OF_NAMES in AutoRandomF*.cc
				 {EQ_NAMES, ISO_NAMES, DOMAIN_NAMES,
				  TYPE_NAMES, MONOTONE_NAMES, MODE_NAMES,
				  OUTPUTMODE_NAMES, REPORTCOORD_NAMES,
				  UNITS_NAMES, COORD_SYS_NAMES,
				  CARTESIAN_SYS_NAMES, TYPEOF_PARAM_NAMES};

//  int SYS_TO_ISO[nr_coord_sys_proj] =  {GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ};

char STANDARDPARAM[MAXPARAM][MAXCHAR],
  STANDARDSUB[MAXSUB][MAXCHAR];

model **KEY() { return KEYT()->KEY; }

int currentRegister() {
  KEY_type *p = KEYT();
  if (p == NULL) BUG;
  return p->currentRegister;
}

void set_currentRegister(int cR) {
  KEY_type *p = KEYT();
  if (p == NULL) BUG;
  p->currentRegister = cR; 
}

//int zaehler = 0;
double *ZERO(int dim, KEY_type *KT) {
  //   printf("zero enter %lu %lu\n", KT, KT->zerox);
  assert(KT != NULL);
  if (dim > KT->nzero) {
    FREE(KT->zerox); 
    KT->nzero = dim;
    KT->zerox = (double*) CALLOC(KT->nzero, sizeof(double));
  }
  //  printf("%lu %lu %d %10g\n", KT, KT->zerox, KT->nzero,  KT->zerox[0]);
  // if (zaehler++) crash();
  return KT->zerox;
}

double *ZERO(model *cov) {
  assert(cov != NULL && cov->base != NULL);
  return ZERO(PREVTOTALXDIM, cov->base); // i_ncol
}

void AtZero(int *info, model *cov, double *v) {
  double *zero = ZERO(cov);
  if (isKernel(PREV)) NONSTATCOV(zero, zero, info, cov, v)
    else COV(zero, info, cov, v);
}

void AtZero(model *cov, double *v) {
  DEFAULT_INFO(info);
  AtZero(info, cov, v);
}
 

bool CheckListmodel(){
  assert((bool) 5 == 1);
  assert(MODEL_MAX == 30); // otherwise change REG_NAMES
  assert(OtherType == 20 && LASTTYPE == 30); // otherwise change TYPE_NAMES, 
  //                                           CAT_TYPE_NAMES[OtherType + 1]
  //                                           Type-Consistency in question.cc
 assert(LAST_ISO == 20); // otherwise change ISO_NAMES
  //  assert(MAXMPPDIM <= MAXSIMUDIM); // ZERO
  // assert(CUTOFF_THEOR == 4);/* do not change this value as used in RFmethods.Rd */

#define ntypefcts 28
  const char *typefcts[ntypefcts] =
    {SYMBOL_PLUS, SYMBOL_MULT, SYMBOL_S, "$power", "fixcov",
     "identity", "M", "matern", "tbm", RM_USER,
     "whittle", "nugget", "null", "shape", "trafo",
     "setparam", "p", "c", "Simulate", "direct", "nuggetIntern",
     "binaryprocess", "gauss.process", "Cov", "CovMatrix",
     "Variogram", "Pseudovariogram", "loglikelihood"};

  int nr = 0,
    err=NOERROR;
  defn *C = NULL;
  for ( ; nr<currentNrCov; nr++) {
    //    printf("n=%d(%d) %s\n", nr, currentNrCov, C->name);

    
    C = DefList + nr; // nicht gatternr    
      
    err = 1;
    if (isManifold(SYSTYPE(C->systems[0], 0)) && C->TypeFct == NULL) 
      goto ErrorHandling;

    err = 2;
    if (isMathDef(DefList + nr)) {
      if (DefList[nr].variants > 2 && nr != CONST) goto ErrorHandling;
    }

    err = 3;
    if ((equalsParamDepI(C->systems[0][0].iso) ||
	 equalsParamDepD(C->systems[0][0].dom))
	&& (C->TypeFct == NULL || C->setDI == NULL)) goto ErrorHandling;

    err = 4;
    if (C->setDI == NULL) {
      err = 104;
       if (isParamDepI(C) || isParamDepD(C)) goto ErrorHandling;
    } else {
      if (!(isParamDepI(C) || isParamDepD(C) || isSubModelI(C)))
	goto ErrorHandling;
    }

    err = 5;
    if (C->TypeFct != NULL) {
      Types type = SYSTYPE(C->systems[0], 0);
      if(//isProcess(type) ||
	 isPointShape(type) || equalsRandom(type)
	 //	 || equalsInterface(type)
	 ) goto ErrorHandling;
    }

    err = 6;
    for (int k=0; k<C->kappas; k++) {
      if (C->kappanames[k][0] == ONEARGUMENT_NAME
	  && C->kappanames[k][1] >= '0'
	  && C->kappanames[k][1] <= '9') goto ErrorHandling;
    }

    err = 7;
    if ( !(!equalsUnreduced(ISO(C->systems[0], 0)) || C->variants == 1 ||
	   (C->variants == 2 && equalsUnreduced(ISO(C->systems[1], 0))) ||
	   ((C->variants >= 3) &&
	    C->Iallowed != NULL && C->Iallowed != allowedIstandard &&
	    C->Iallowed != allowedPrevModelI &&
	    equalsIsotropic(ISO(C->systems[1], 0)) &&
	    equalsEarthIsotropic(ISO(C->systems[2], 0)) &&
	    (C->variants == 3 || equalsUnreduced(ISO(C->systems[3], 0))))
	   )) goto ErrorHandling;

    err = 8;
    if (C->TypeFct != NULL && C->name[0] != '-') {
      int i;
      for (i=0; i<ntypefcts; i++) if (!STRCMP(typefcts[i], C->name)) break;
      if (i >=ntypefcts) {
	PRINTF("Martin, read instructions for typefcts: '%s'\n", C->name);
	/* 
	   instructions: note that is** in questions.h has a different
	   behavious if TypeFct is given. Check and enter in the list if ok.

	   27.4.2019: isDefCL has been rewritten recently so that this point
	   may not exist anymore
	 */
	//	goto ErrorHandling;
      }
    }

    err = 9;
    if (nr != COVARIATE && nr != FIXCOV && nr != VARIOGRAM2COV &&
	nr != VAR2COV_PROC && nr != PREDICT_CALL) {
      //      printf("%s\n", C->name);
      for (int k=0; k<C->kappas; k++) {
	if (C->kappatype[k] == VECSXP) {
	  PRINTF("Martin, see KeyInfo.cc fctn Param to deal with %s\n",C->name);
	  goto ErrorHandling;
	}
      }
    }

    err = 10; // MISMATCH=-4
    // 
    Types type = SYSTYPE(C->systems[0], 0);
    if (nr > LAST_TRAFO && C->maxmoments < 0 && C->maxmoments != SUBMODEL_DEP &&
	C->maxsub != 0 && C->Init != init_failed && !isInterface(type)) {
      if (C->cov != ave && C->cov != Id && C->check != checkstp  //&&
	  // C->check != check_mcmc_pgs && C->check != check_Zhou
	  ) {
	PRINTF("%s %d nr=%d\n", C->nick, C->maxmoments, nr);
	goto ErrorHandling;
      }
    }
 
    err = 11;
    if (nr != VARIOGRAM2COV && nr !=  VAR2COV_PROC && nr != PREDICT_CALL) {
      for (int k=0; k<C->kappas; k++) {
	if (C->kappatype[k] == VECSXP &&
	    (k==C->kappas-1 || STRCMP(C->kappanames[k+1], COVARIATE_RAW_NAME))){
	  PRINTF("%s: Martin, all functions that have additional points inside must set KT->rawConcerns = neverRaw; then included in the exception list", C->name);
	  goto ErrorHandling;
	}
      }
    }
  }

  printf("urgent to do: replace CORES by KEYT()->global_utils \n");// OK
  printf("partially match of location with kriging\n"); // OK
  printf("!!!! ACHTUNG ! bei Verwendung von TALLOC auf Sextra basierend, darf bis zu END_TALLOC kein cov stehen, oder es muss sichergestellt werden  dass die Fkt mit argument 'cov' nicht auch noch TALLOC aufruft!\n"); // OK
  printf("fft is currently disabled\n"); //
  printf("smith model \n"); // OK
  printf("steht was wichtiges in famillies.cc.spaeter?\n"); // OK
  printf("eps<1.0 in huetchen.cc!"); // OK
  printf("rpbernoulli auch fuer prozesse"); // OK
  printf("brownresnick ueberarbeiten. insb. warum dort TransformLoc aufgerufen wrid"); // OK
  printf("trafo in variogram.cc nach vorne ziehn, so dass Gitter bleiben wenn nur gestreckt, also statt grid_expand=True ein AVOID"); // OK
  printf("case !grid, proj<0, separable to be programmed"); // OK
  printf(" Doku rfsimulateadvanced "); // OK
  printf(" ength(err.model) > 1 in RFinterpolate R; xRMranef in RMmodelsSpecial.R und kriging.R"); // OK
  printf("RFgui "); // OK
  printf(" RMS: quasi projektmatrizen: blockmatrizen"); // OK
 printf(" sequential: 'back' in Sequential fuer whittle (bei nu=markov field)exakt suchen und als Specific implementieren? -- pref umsetzen in whittle etc"); // OK
 printf(" gausslikeli ruft covvario auf und x wird gegebenefalls immer wieder von neuem transformiert"); // OK
 printf(" krigin variance"); // OK
 printf(" USE_OWN_ALG + GAUSS * Chol -> Ux + test in utils.cc ausprobieren; http://www.netlib.org/blas/"); // OK
 printf(" "); // OK
  //  printf("done\n");
  return true;


 ErrorHandling:
  // if (C != NULL) PRINTF("%s: ", C->nick);
  PRINTF("Fehler in CheckListmodel: err=%d nr=%d (%s) %d\n", err, nr,
	 DefList[nr].nick, C == NULL);
  
  // PRINTF("\n\nFehler in CheckListmodel: name = %s %d (err=%d)\n\n",
  //	 C == NULL ? "xxx" : C->nick, nr, err);
  return false;
}

void InitModelList() {
  assert(currentNrCov == UNSET); // otherwise something went wrong with the call
 
  for (int i=0; i<MAXPARAM; i++) SPRINTF(STANDARDPARAM[i], "k%d", i+1);
  for (int i=0; i<MAXSUB; i++) SPRINTF(STANDARDSUB[i], "u%d", i+1);
  /* ja nicht setzen !! macht riesen aerger, da RF opt ions InitModel
    nicht aufruft:
    for (i=0; i<MAX UNITS; i++) {
    STRCPY(glo bal->general.newunits[i], "");
    STRCPY(glo bal->general.curunits[i], "");
    STRCPY(glo bal->general.varunits[i], "");
  }
  */


   // init models

  if (DefList!=NULL) {
    PRINTF("List of covariance functions looks already initiated.\n"); 
    return;
  }
  DefList = (defn*) MALLOC(sizeof(defn) * (MAXNRCOVFCTS+1));
  // + 1 is necessary because of COVINFO_NULL that uses the last + 
  currentNrCov = 0;

 
  // *******************
  // **** RO-models ****
  // *******************

   
  FIRSTGATTER0 =  // 0 -- always first
    IncludeModel("#",  OtherType, 1, 1, 0, NULL, PREVMODEL_D, PREVMODEL_I,
  		 checkNotOK, NULL, PREF_NOTHING, true, SUBMODEL_DEP,
  		 SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
 
 assert(FIRSTGATTER0 == FIRSTGATTER);
  addCov(stat2, D_2, DD_2, D3_2, D4_2, inverse2, inverse_nonstat2);
  addCov(nonstat2);// 
  addlogCov(logstat2, nonstat_log2, inverse_log_nonstat2);
  // addCov(iso2iso, D_2, DD_2, inverse2, inversenonstat2);
  RandomShape(INFTY, struct2, init2, do2, dorandom2, true, true, false); 
  
  ISO2ISO0 = addFurtherCov(ErrCov, ErrD, ErrD); // 1
  SP2SP0 = addFurtherCov(ErrCov, ErrD, ErrD); // 2
  SP2ISO0 = addFurtherCov(ErrCov, ErrD, ErrD); // 3
  S2ISO0 = addFurtherCov(ErrCov, ErrD, ErrD); // 4
  S2S0 = addFurtherCov(ErrCov, ErrD, ErrD);// 5
  SId0 = addFurtherCov(ErrCov, ErrD, ErrD);// 6
  S2SP0 = addFurtherCov(ErrCov, ErrD, ErrD);// 7
  E2EIso0 = addFurtherCov(ErrCov, ErrD);// 8
  E2E0 = addFurtherCov(ErrCov, ErrD);// 9
  E2SphIso0 = addFurtherCov(ErrCov, ErrD);// 10
  E2Sph0 = addFurtherCov(ErrCov, ErrD);// 11
  Sph2SphIso0 = addFurtherCov(ErrCov, ErrD);// 12
  Sph2Sph0 = addFurtherCov(ErrCov, ErrD);// 13

  LASTGATTER0 = Sph2Sph0;
  assert(LASTGATTER0 == LASTGATTER);
  
  FIRST_TRAFO = EARTHKM2CART= // 14
      IncludeModel(">",  OtherType, 1, 1, 0, NULL,
		   PREVMODEL_D, PREVMODEL_I, // dummy values
		   checkEarth, NULL, PREF_NOTHING, true, SUBMODEL_DEP,
		   4, (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  
  addCov(EarthKM2CartStat);
  addlogCov(EarthKM2Cart);//
 
  EARTHMILES2CART = addFurtherCov(EarthMiles2CartStat, ErrD);// 15
  addlogCov(EarthMiles2Cart);// 

  FIRST_PLANE = EARTHKM2GNOMONIC =
    addFurtherCov(Earth2GnomonicStat, ErrD);// 16
  addlogCov(Earth2Gnomonic);// 

  EARTHMILES2GNOMONIC =  CopyModel(">", EARTHKM2GNOMONIC); // 17

  EARTHKM2ORTHOGRAPHIC = addFurtherCov(EarthKM2OrthogStat, ErrD);// 18
  addlogCov(EarthKM2Orthog);// 

  EARTHMILES2ORTHOGRAPHIC = addFurtherCov(EarthMiles2OrthogStat, ErrD);// 19
  addlogCov(EarthMiles2Orthog);// 
 
  LAST_PLANE = EARTHMILES2ORTHOGRAPHIC;
  LAST_TRAFO =  EARTHMILES2ORTHOGRAPHIC;
  assert(LAST_TRAFO == 19);

  pref_type pplus =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 0, 5};
  //                  CE CO CI TBM Sp di sq Tr av n mpp Hy spf any
  PLUS = 
    IncludeModel(SYMBOL_PLUS, ManifoldType, 1, MAXSUB, 1, NULL,
		 PREV_SUB_D, PREV_SUB_I,
		 checkplus, rangeplus, pplus,  false,
		 SUBMODEL_DEP, SUBMODEL_DEP,
		 (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  // Achtung in covmatrix_plus wird SELECT_SUBNR verwendet!
  nickname(PLUS_NICK);
  kappanames("trend", INTSXP);
  addCov(plus, Dplus, DDplus);
  addCov(nonstat_plus);
  addTBM(NULL, spectralplus);
  RandomShape(0, structplus, initplus, doplus);
  addReturns(covmatrix_plus, iscovmatrix_plus);
  setptwise(pt_submodeldep);
  addTypeFct(Typeplus);
  setDI(allowedDplus, allowedIplus, NULL);
  
  
  pref_type pmal =  {5, 0, 0,  5, 0, 5, 5, 0, 0, 0, 0, 0, 4, 5};
  //                 CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  MULT =
    IncludeModel(SYMBOL_MULT, ManifoldType, 1, MAXSUB, 0, NULL, PREV_SUB_D, PREV_SUB_I,
		 checkmal, NULL, pmal, false, SUBMODEL_DEP, SUBMODEL_DEP,
		 (ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  nickname(MULT_NICK);
  addCov(mal, Dmal);
  addCov(nonstat_mal);
  addlogCov(logmal, nonstat_logmal, NULL);
  RandomShape(0, structmal, initmal, domal);
  setptwise(pt_submodeldep);
  addTypeFct(Typemal);
  setDI(allowedDplus, allowedIplus, NULL);


  pref_type pS=  {5, 0, 0,  5, 5, 5, 5, 0, 0, 5, 0, 0, 1, 5};
  //              CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  FIRSTDOLLAR = IncludeModel(SYMBOL_S,  ManifoldType, // to do: tcftype durch einen allgemeinen Type ersetzen, da auch Trend dem "$" folgen kann. Z.Z. nicht moeglich.
			1, 1, 5, kappaS, // kappadollar,
			     PREV_SUB_D, PREV_SUB_I, checkS, rangeS, pS,
			false, SUBMODEL_DEP, SUBMODEL_DEP,
			(ext_bool) SUBMODEL_DEP, MON_SUB_DEP);
  // do not change Order!!
  nickname(S_NICK);
  kappanames("var", REALSXP, "scale", REALSXP, "anisoT", REALSXP,
	     "Aniso", REALSXP, "proj", INTSXP);
  change_typeof(DVAR, RandomOrShapeType);
  change_typeof(DSCALE, RandomOrShapeType);
  change_typeof(DAUSER, ShapeType);
  subnames("phi");
  addCov(Siso, DS, DDS, D3S, D4S, inverseS, inversenonstatS); // unterscheidung nur wegen der 
  //  geschwindigkeit, also Siso ist ein sehr haeufiger Spezialfall von Sstat
  addCov(nonstatS);
  addlogCov(logSiso, NULL, inverse_log_nonstatS);
  addLocal(coinitS, ieinitS);  
  addTBM(tbm2S, NULL, spectralS);
  nablahess(nablaS, hessS);
  RandomShape(INFTY, structS, initS, doS, true, true, false);
  addReturns(covmatrixS, iscovmatrixS);
  Taylor(RF_NA, RF_NA, RF_NA, RF_NA);
  TailTaylor(RF_NA, RF_NA, RF_NA, RF_NA);
  setptwise(pt_submodeldep); 
  addTypeFct(TypeS);
  setDI(allowedDS, allowedIS, NULL); //setS
 
   
  LASTDOLLAR = addFurtherCov(Sstat, DS, DDS); // 20.8.14 aus ErrCov (wieder)
  //                                        D2 gemacht
  addCov(nonstatS);
  addlogCov(logSstat, nonstat_logS, NULL);

  // printf("%d\n",  currentNrCov); BUG;

  pref_type pPowS=  {5, 0, 0,  5, 5, 5, 5, 0, 0, 5, 0, 0, 1, 5};
  //                CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  POWER_DOLLAR = 
    IncludeModel("$power", ManifoldType, // to do: tcftype durch einen allgemeinen Type ersetzen, da auch Trend dem "$" folgen kann. Z.Z. nicht moeglich.
		 1, 1, 3, NULL, // kappadollar,
		 PREV_SUB_D, PREV_SUB_I, checkPowS, rangePowS, pPowS,
		 true, SUBMODEL_DEP, SUBMODEL_DEP, (ext_bool) SUBMODEL_DEP,
		 MON_SUB_DEP);
  // do not change Order!!
  nickname("Spower");
  kappanames("var", REALSXP, "scale", REALSXP, "pow", REALSXP);
  subnames("phi");
  addCov(PowS, NULL, NULL, inversePowS, inversenonstatPowS); // unterscheidung nur wegen der 
  //  geschwindigkeit, also Siso ist ein sehr haeufiger Spezialfall von Sstat
  addCov(nonstatPowS);
  addlogCov(logPowS, nonstat_logPowS, NULL);
  // addLocal(coinitS, ieinitS);  
  RandomShape(INFTY, structPowS, initPowS, doPowS, true, true, true);
  Taylor(RF_NA, RF_NA, RF_NA, RF_NA);
  TailTaylor(RF_NA, RF_NA, RF_NA, RF_NA);
  addTypeFct(TypePowS);
  setptwise(pt_submodeldep); 
  

  ////////////////////////////////////////////////////////////
  includeCovModels();
  //  includeAsymmetricModels();
  includeOtherModels();
  ////////////////////////////////////////////////////////////

   IncludeModel("minus", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMath, PREF_TREND,
	       false,SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP, "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathminus);
  AddVariant(TrendType, PREVMODEL_I);
  set_type(DefList[currentNrCov-1].systems[0], 0, ShapeType);
 
  IncludeModel("plus", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMath, PREF_MATHDEF, 
	      false,SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP, "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathplus);
  AddVariant(TrendType, PREVMODEL_I);
  setptwise(pt_submodeldep); 
 
  MATHDIV = IncludeModel("div", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkDivMult, rangeMath, PREF_MATHDEF, 
	      false,SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP,  "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathdiv);
  AddVariant(TrendType, PREVMODEL_I);
  setptwise(pt_submodeldep); 
 

  IncludeModel("mult", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkDivMult, rangeMath, PREF_MATHDEF, 
	      false,SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "y", REALSXP,  "factor", REALSXP);
  change_sortof(MATH_FACTOR, TRENDPARAM);
  addCov(Mathmult);
  setptwise(pt_submodeldep); 
  AddVariant(TrendType, PREVMODEL_I);

  CONST =
    IncludeModelR(R_CONST, MathDefType, 0, 0, 2, NULL, XONLY, PREVMODEL_I,
		 check_c, rangec, PREF_MATHDEF, 
		 false, SCALAR, PREVMODEL_DEP, falsch, NOT_MONOTONE);
  kappanames(CONST_A_NAME, REALSXP, COVARIATE_NAME_NAME, STRSXP);
  change_sortof(CONST_C, TRENDPARAM);
  change_sortof(CONST_NAME, IGNOREPARAM);
  addCov(Mathc);
  AddVariant(TrendType, PREVMODEL_I);
  AddVariant(NegDefType, PREVMODEL_I);
  AddVariant(TcfType, PREVMODEL_I);
  setDI(NULL, allowedItrue, NULL);
  setptwise(pt_optionsdep); 

  PROJ_MODEL = 
  IncludeModel(R_P, MathDefType, 0, 0, 4, NULL, XONLY, PARAMDEP_I,
	       checkproj, rangeproj, PREF_MATHDEF, 
	       INTERN_SHOW,  SCALAR, INFDIM-1, falsch, NOT_MONOTONE);
  kappanames("proj",  INTSXP, "new", INTSXP, "factor",  REALSXP,
	     COVARIATE_NAME_NAME, STRSXP);
  change_typeof(PROJ_ISO, NN2);// "new"
  change_sortof(PROJ_FACTOR , TRENDPARAM);
  change_sortof(PROJ_NAME, IGNOREPARAM); 
  addCov(proj);
  AddVariant(TrendType, PREVMODEL_I);
  setDI(NULL, allowedIp, setproj);
  addTypeFct(Typeproj);
  
  BIND = 
  IncludeModel(R_C, MathDefType, 0, 0, 18, NULL, PREV_SUB_D, PREV_SUB_I,
	       check_bind, rangeMath, PREF_TREND, 
	       INTERN_SHOW, PARAM_DEP, 1, falsch, NOT_MONOTONE);
  kappanames("a", REALSXP, "b", REALSXP, "c", REALSXP, 
	     "d", REALSXP, "e", REALSXP, "f", REALSXP,
	     "g", REALSXP, "h", REALSXP, "i", REALSXP,
	     "j", REALSXP, "l", REALSXP, "m", REALSXP,
 	     "n", REALSXP, "o", REALSXP, "p", REALSXP, 
 	     "q", REALSXP, // R: 16, C: 15
	     "ncol", INTSXP, "factor", REALSXP);
  change_sortof(DefList[BIND].kappas - 1, TRENDPARAM);
  addCov(Mathbind);
  AddVariant(TrendType, SUBMODEL_I);
  assert(BIND_VARIABLES == 16);
  assert(DefList[BIND].kappas == BIND_VARIABLES + 2);
  set_type(DefList[currentNrCov-1].systems[0], 0, ShapeType);
  setDI(allowedDbind, allowedIbind, NULL);
  addTypeFct(Typebind);

 
  IncludeModel("is", MathDefType, 0, 0, 3, NULL, XONLY, PREVMODEL_I,
	       checkMath, rangeMathIs, PREF_TREND, 
	      false, SCALAR, 1, falsch, NOT_MONOTONE);
  kappanames("x", REALSXP, "is", INTSXP, "y", REALSXP);
  change_typeof(IS_IS, NN1);
  addCov(MathIs);
  AddVariant(TrendType, PREVMODEL_I);
  set_type(DefList[currentNrCov-1].systems[0], 0, ShapeType);
  setptwise(pt_posdef); 

  includeStandardMath();
  
  assert(CheckListmodel());

}
