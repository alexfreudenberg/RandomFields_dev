

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

#ifndef RFsimu_H
#define RFsimu_H 1


// ACHTUNG : REIHENFOLGE WICHTIG. def.h zuerst; obige defines immer oben
#include "def.h"
#include "init.h"
#include <inttypes.h>
typedef unsigned int Uint;
typedef uint64_t  Ulong;
typedef int64_t  Long;


// ohne: 0.056
// mit: 0.268
// ohne ASSERT_GATTER : 0.199
// mit dummy ASSERT_GATTER : 0.194


#define showfree !true 
#define DOPRINT true

// in /home/schlather/TMP/RandomFieldsUtils/include/utils.h:
// #define MEMCOPY(A,B,C) memory_copy(A, B, C)




// covariate_storage *S = cov->Scovariate;
// M-x replace-regexp RET c[ad]+r RET \&-safe RET
// grep "\(.*\)_storage.*\*\(.*\)=.*->S.*" getStorage(\2, \1); 



#include <string.h>
#include <General_utils.h>
#include "AutoRandomFields.h"
#include "auxiliary.h"
#include "RandomFields.h"
#include "xport_import.h"
#include "PoissonPolygon.h"
#include "Options.h"

//#undef assert
//#define assert(X) // runter auf 0.115 


#ifndef RFERROR
#define RFERROR error
#endif

// never use Rpri ntf/PRINTF/print within DO_PARALLEL !!!!!!
#undef print //
#define print NEVER_USE_print_or_PRINTF_WITHIN_PARALLEL

#define MAX_LEN_EXAMPLES 4
#define EIGENVALUE_EPS 1e-15  // used in GetTrueDim
#define GENERAL_PRECISION 5e-15

///////////////////////////////////////////////////////////////////////
// BASIC DIMENSIONS AND VARIABLES
///////////////////////////////////////////////////////////////////////
#define PARAMMAXCHAR MAXCHAR
#define MAXLONGCHAR 40 // for RFoptions that are string

//#define MAX_MLE_ELMNTS 10
/* max number of submodels -- NEVER MORE THAN 255 */
#define MAXTBMDIM 3
#define MAXTAYLOR 3
#/* if explicite simulation method is not given by ????????
   user, the programme will first try a methods 
   according to "first"; 
   then it tries the other methods according to
   a default list: (CircEmbed only if the simulation is on a grid)
   1 dim: CircEmbed, Direct
   2 dim: CircEmbed, TBM, spectralTBM, 
   AdditiveMpp, Direct
   3 dim: CircEmbed, TBM, Direct
   the method actually used is stored in actually_used. 
*/ 
#define MAXLILIGRIDDIM 10


#define NATSCALE_EXACT 1 
#define NATSCALE_ORNUMERIC 2
//#define NATSCALE_APPROX 2
#define NATSCALE_MLE 3 /* check fitvario when changing; best NOT to change !! */

#define LISTOF 1000  /* may NOT interfere with define of *SXP in Rinternal.h*/

#define ONEARGUMENT_NAME 'k'
#define DEFAULT_SUBNAME 'C'

   

///////////////////////////////////////////////////////////////////////
// Preferences for the choice of the gaussian simulation method
#define PREF_BEST 5
#define PREF_NONE 0
#define PREF_NUGGET 1
#define PREF_PENALTY Nothing
#define PREF_FACTOR Nothing
// larger than Pref * Nothing
#define LOC_PREF_NONE -1000



///////////////////////////////////////////////////////////////////////
#define I_COL_NA -1


///////////////////////////////////////////////////////////////////////
// Kappas
#define SIZE_NOT_DETERMINED 0 // don't change the value !
#define OUT_OF_RANGE -1


//////////////////////////////////////////////////////////////////////
// the different levels of printing

#ifndef PL_IMPORTANT // SEE utils.h in RandomFieldsUtils
#define PL_IMPORTANT 1 // crucial messages
#define PL_BRANCHING 2 // user relevant RPgauss branching etc
#define PL_DETAILSUSER 3
#define PL_RECURSIVE 4 
#define PL_STRUCTURE 5 // see also initNerror.ERROROUTOFMETHOD
#define PL_ERRORS  6 // only those that are caught internally

#define PL_FCTN_DETAILS 7  // R
#define PL_FCTN_SUBDETAILS 8

#define PL_COV_STRUCTURE 7 // C
#define PL_DIRECT_SEQU 8
#define PL_DETAILS 9
#define PL_SUBDETAILS 10
#endif



///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////
// type of covariance functions that need distinguished treatment
// Reihenfolge beachten! Hoehere Nummern sind meist in kleinere umwandelbar


// vdim is "strictly" passed from bottom to top
#define SCALAR 1 

// way of implementing simulation methods:
#define NOT_IMPLEMENTED 0 /* don't change this val_ue except to false */
#define IMPLEMENTED 1     /* don't change this value except to true */
#define NUM_APPROX 2 
#define GIVEN_METH_IGNORED 3
#define HYPERIMPLEMENTED 4



///////////////////////////////////////////////////////////////////////
// BASICS
///////////////////////////////////////////////////////////////////////

 
typedef enum usekernel_type {noKernel, allowKernel,forceKernel} usekernel_type;

typedef enum ext_bool {falsch = false, 
		       wahr = true, 
		       Paramdep=PARAM_DEP,  // necessarily negative !!
		       Submodeldep=SUBMODEL_DEP, // necessarily negative !!
		       Huetchenownsize=2,
		       BothOK = 3
                     } ext_bool;

typedef enum TaylorCoeff {TaylorConst, TaylorPow, TaylorExpConst, 
			  TaylorExpPow} TaylorCoeff;
// Const * x^Pow * e^{-ExpConst * x^ExpPow} //
// Falls finite range, so wird in TaylorConst der range abgelegt



typedef double **coord_type; 
struct location_type {
  int 
    timespacedim,      /* total dimension, including time and space */
  //length[M A XS IMUDIM], /* grid==true: number of pixels of each direction 
  //		       grid==false: length[0]=number of spatial points in total
  //		       	     */
    spatialdim, 
    xdimOZ, // physical xdim, without time
    len; // nur beim ersten Element sinnvoll. Aber der 
  // Einfachheit nicht nochmals einen strukt drumrumgemacht, da sowieso
  // meist == 1
 
  int // deleted! lx, lY, -- 14.1.21
  // NEVER USE THIS EXCEPT IN RESETTING BY loc_set
  /*          
	      the total number of locations of the user's x vector, i.e.
	      * except for grids where lx=3
              * physical length = xdimOZ * lx * (lx - 1) /2 for distances
	      * physical length = xdimOZ * lx else
	     */
    spatialtotalpoints,
    totalpoints,     /* 14.1.21 uebernimmt die Aufgabe von vorher lx,
			um anzuzeigen, dass nichts gesetzt ist.
			number of pixels in total;
			
			    == lengthx * lengthy * lengthz if grid==true
			 */
    spatialtotalpointsY,
    totalpointsY,  /* number of pixels in total for conditioning points */
    rawset, *rawidx;  /* for subset of the locations, e.g. from fitgauss
			 or neighbourhood kriging;
			 set of locations may have split into various 
			 subsets. So i_idx and set do not refer to the 
			 correct indices for data in 'covariate' and 'fix'.
			 The correct indices are given here by transformations
		      */
  bool grid, gridY, /* simulation on a grid required by user? */
    delete_x, delete_y,  /* is loc->x only a reference? -- used in MLE */
    distances,//are distances given, and not locations ? -- mostly for matrices
    Time;             /* is a time component given ? */

  coord_type xgr, grY;      /* the coordinates */  

  double *x,   /* sortiert 1. Koord. 1. Punkt, 2. Koord 1. Punkt, ... */
    T[3],      /* gridtriple definition for the time components --
		  doubled in last coordinate of *x
	       */
    *Y, /* only used by Covariance function evaluation */
    TY[3], /* also only used by covariance fctn evaluation */
    *caniso;   /* only used for grid. Cummulates the 
		  anisotropy matrices when initS is called
	       */
 int cani_ncol, cani_nrow;
};




///////////////////////////////////////////////////////////////////////
// GENERAL PARAMETERS FOR THE SIMULAIONMETHODS
///////////////////////////////////////////////////////////////////////

typedef
struct model model;
typedef
struct gen_storage gen_storage;
typedef int (*structfct)(model *cov, model **newcov);
typedef int (*initfct)(model *cov, gen_storage *s);
typedef void (*dofct)(model *cov, gen_storage *s);
typedef void (*finaldofct)(model *cov, double *res, int n, gen_storage *s);
typedef void (*do_random_fct)(model *cov, double *v);
typedef void (*param_set_fct)(model *to, model *from, int variant);
typedef Types (*type_fct)(Types required, model *cov, isotropy_type i);


typedef enum matrix_type {//TypeMid,
                          TypeMiso, TypeMdiag, 
			  TypeMtimesepproj, // TypeMtimesep and TypeMproj
			  TypeMtimesep, // last column is zero, but last entry
			  TypeMproj, // including scale: values need not be 1 
			  TypeMany} matrix_type;

struct simu_storage {
  bool active, pair;   /* has the init_procedure been called successfully? */
  int expected_number_simu;
};



///////////////////////////////////////////////////////////////////////
// COVARIANCE SPECIFICATIONS
///////////////////////////////////////////////////////////////////////


typedef int pref_shorttype[Nothing + 1];
struct simple_range_type {
  double min;
  double max;
  bool openmin;
  bool openmax;
  double pmin;
  double pmax;
};
struct range_type {
  double min[MAXPARAM];
  double max[MAXPARAM];
  bool openmin[MAXPARAM];
  bool openmax[MAXPARAM];
  double pmin[MAXPARAM];
  double pmax[MAXPARAM];
};

typedef double *param_type[MAXPARAM];
struct sexp_type {
  bool Delete;
  SEXP sexp;
};

struct listoftype {
  bool deletelist;
  double **lpx;
  int Rtype, // LISTOF + REALSXP
    len, // identical to nrow of calling structure
    *ncol, // in case of random parameter, they might be given
    *nrow; // although p[.] is NULL; latter set in CHECK
};

typedef double *naptr_type[MAX_NA];
typedef model *covptr_type[MAX_NA];
//typedef int *elptr_type[MAX_MLE_ELMNTS];
// typedef double *internal_type[MAXINTERNALPARAM];

typedef double (*spectral_density)(double *, model *cov);

// TO DO UEBERPRFEN : ALLE *_storage duerfen nur waeren check/init
// beschrieben werden! Ansonsten ist TALLOC, siehe z.b. tbm in operator.gauss.cc
// zu verwenden !!
struct spectral_storage {
  //init:
  spectral_density density;
  double sigma, E[MAXTBMSPDIM];
  int nmetro;
  double sub_var_cum[MAXSUB]; 
  // simu:
  double phistep2d, phi2d, prop_factor;
  bool grid;
  
};


typedef enum user_given_type {ug_explicit, ug_implicit, ug_internal}
  user_given_type;


struct mle_storage {
  naptr_type MEMORY;
  covptr_type COVMODELS;
  double *PT_VARIANCE;
  int NAS;
};



#define MAXLOCALINSTANCES 3
struct localfactstype {
  int instances;
  int msg[MAXLOCALINSTANCES];
  double value[MAXLOCALINSTANCES];
};

typedef bool allowedD_type[LAST_DOMAINUSER + 1];
typedef bool allowedI_type[LAST_ISOUSER + 1];


#define INFO_N_X 0
#define INFO_IDX_X 1
#define INFO_EXTRA_DATA_X 2
#define INFO_N_Y 3
#define INFO_IDX_Y 4
#define INFO_EXTRA_DATA_Y 5 
#define INFOS_PER_COORD 3
#define INFO_SIZE (2 * INFOS_PER_COORD)
#define INFO int VARIABLE_IS_NOT_USED *info
#define DEFAULT_INFO(info)						\
  int info[INFO_SIZE] = {1, NA_INTEGER, false, 1, NA_INTEGER, false}


typedef void (*rangefct)(model *cov, range_type* ra);
typedef void (*rangefct_multi)(model *cov, int, int, int,
			       simple_range_type* ra);
typedef int (*checkfct)(model *cov); 
typedef void (*covfct)(double *, int *, model*, double*);/* h, cov, result*/
typedef void (*nonstat_covfct)(double *, double*, int *,
			      model*,  double*); /* x,y, cov, result */
typedef void (*logfct)(double *, int *, model*, double* v, double* Sign);
typedef void (*nonstat_logfct)(double *, double*, int *,
			      model*,  double* v, double* Sign); 
typedef void (*inverse_fct)(double *v, model *cov, double *x); /* parameters,; natural scaling */
typedef void (*nonstat_inv)(double *, model*, double*, 
			       double*); /* fctvalue cov, inv.x inv.y */
//typedef void (*aux_covfct)(double *, double*, int*, double, // obsolete ?
//			      model*,  double*); /* x,y, Aux, cov, result */
typedef void (*return_fct)(model*, double*); /* cov, result */ 
typedef void (*return_fctalpha)(double, model*, double*); /* cov, result */ 
typedef void (*tworeturns_fct)(model*, double*, double*); /* cov, result */ 
typedef void (*return_covmat)(model*, bool, double*); /* cov, result */ 
typedef char (*ext_bool_ret_fct)(model*); /* covt */ 
typedef void (*getlocalparam)(model *, localfactstype *);
typedef bool (*altlocalparam)(model *);
typedef void (*minmaxfct)(model *, double *); 


typedef void (*spectral_do)(model *, gen_storage *, double *);

typedef void (*draw_random) (model *cov, double *random);
typedef double (*log_mixdens)(double *x, double logV, model *cov);

typedef void (*sd_fct)(gen_storage *s, model *cov);
		       				   
typedef int (*hyper_pp_fct)(double, double*, double*, model *, bool, 
			    double**, double**, double**);
typedef void (*size_fct)(int i, model *cov, int *nr, int *nc);
typedef sortsofparam  (*sortof_fct)(model *cov, int k, int row, int col,
				    sort_origin origin);

typedef bool (*setDI_fct)(model*); /* false, falls auf PREVISO angewiesen
				      und nicht gesetzt 
				   */ 
typedef bool (*allowedD_fct)(model *cov); /* function returns false if a 
					     selection of models is allowed;
					     the selection is given in 
					     cov->Dallowed. It return true 
					     if none the options should not be
					     tested one after each other, but
					     PREVDOM should be instead.
					   */
typedef bool (*allowedI_fct)(model *cov);



struct system_type {
  int nr,  /* cov number of trafo bzw. eigene Nummer bei own */
    last,
    logicaldim, 
    maxdim, /* maxdim of model combined with information of the submodels
	       this can be UNSET for certain models
	     */                 
    xdim,  /* current xdim *including time* (when entering # or >),
	      which may differ from logicaldim for isotropic or
	      space-isotropic models   */
    cumxdim /* 
	      cummulative cumxdim, including the current position
	     */
     ;
  Types type;  /* DefList->Type may change in case of operators */
  domain_type dom;
  isotropy_type iso;
};

typedef system_type Systems_type[MAXSYSTEMS];

struct defn {  
  char name[MAXCHAR], nick[MAXCHAR],
    kappanames[MAXPARAM][PARAMMAXCHAR],
    subnames[MAXSUB][PARAMMAXCHAR];
  int kappas, // number of parameters  
    minsub, maxsub, variants,
    vdim, maxdim, maxmoments, 
    implemented[Forbidden], Specific, internal, F_derivs, RS_derivs;
  monotone_type Monotone;
  bool primitive,
    subintern[MAXSUB]; // do any subnames match exactly a parameter name?
  ext_bool finiterange;

  Systems_type systems[MAXVARIANTS];
  type_fct TypeFct;
  allowedD_fct Dallowed;
  allowedI_fct Iallowed;  

  SEXPTYPE kappatype[MAXPARAM]; // REALSXP, VECSXP, etc, what is expected
  //                                  from user
  Types kappaParamType[MAXPARAM];//
  //                      RandomType : parameter might be random
  //                      ShapeType : parameter must be fixed
  //                      NN2 : special purpose type
  const char **kappaParamTypeNames[MAXPARAM];
  size_fct kappasize;  // function that gives the size of the argument;
  //                     matrix or vector or scale 
  //                     (all represented as respecitve matrices)
  //
  sortsofparam sortof_tab[MAXPARAM]; // the constants used by the
  //                                     standard paramtype
  sortof_fct sortof; // particularly for MLE, function returning
  //                                   VARPARAM : ANYPARAM,
  //                          but also FORBIDDENPARAM, 

  rangefct range;
  rangefct_multi range_multi;
  checkfct check;
  ptwise_type ptwise_definite;
  pref_shorttype pref;
 
  covfct cov, D, D2, D3, D4, tbm2, nabla, hess, random, logD; //Vtlgen 
  logfct log;
  inverse_fct inverse;
  nonstat_covfct nonstat_cov, nonstat_D, nonstat_random;
  nonstat_inv inverse_nonstat, nonstat_loginverse, inverse_nonstat_D;
  nonstat_logfct nonstatlog;
  param_set_fct param_set;
  //  aux_covfct aux_cov; // complicated cov-model that can be used only
  //  //                     as submodels and that needs an auxiliary argument
  //  //                     for the evaluation
  getlocalparam coinit, ieinit; // set within primitives
  altlocalparam alternative; //getparam: guess for good local 
    // param (cutoff, intrinsic); alternative gives alternative in a 
    // second or third try (used by co and Stein)

  spectral_do spectral;

  draw_random drawmix;
  log_mixdens logmixdens;
 
  structfct Struct;
  initfct Init;
  dofct Do;
  finaldofct FinalDo;
  do_random_fct DoRandom;
  minmaxfct minmaxeigenvalue;
 
  hyper_pp_fct hyperplane;      // hyperplane tessellations         

  
  double Taylor[MAXTAYLOR][TaylorPow + 1], 
    Tail[MAXTAYLOR][TaylorExpPow + 1]; 
  int TaylorN, TailN;

  return_covmat covmatrix;
  //  tworeturns_fct inversecovmatrix;
  ext_bool_ret_fct is_covmatrix;
  /* return value 0 : undefined
                  1 : defined
                >=2 : defined & StandardCall will be very disadvantageous
  */

  setDI_fct setDI;
};







	 


///////////////////////////////////////////////////////////////////////
// STORAGES FOR SPECIFIC MODELS
///////////////////////////////////////////////////////////////////////

// see circembed.cc
struct FFT_storage {
  double* work;
  int *iwork, nseg, maxf[MAXCEDIM], kt[MAXCEDIM], m_fac[MAXCEDIM],
    NFAC[MAXCEDIM][21];
};

struct ce_storage {
  int m[MAXCEDIM], trials,
    halfm[MAXCEDIM], nn[MAXCEDIM], cumm[MAXCEDIM+1], 
    cur_square[MAXCEDIM], max_squares[MAXCEDIM], /* !!!! **** */ 
    vdim; //  added by PM 12/08

  Long mtot, square_seg[MAXCEDIM];
  double **c, **d, smallestRe, largestAbsIm, *aniso;
  complex *gauss1, *gauss2;
  bool positivedefinite,
    stop,
    new_simulation_next,
    cur_call_odd,
     dependent; // eigentlich braucht es nicht waehrend der initialisierung
    // festgelegt zu werden. Ist aber wesentlich einfacher zu handhaben,
    // da sonst bei internal_dosimulate die parameter alle RFparameter alle
    // nochmals gesetzt werden muessen
#ifdef CE_PARALLEL
  FFT_storage *FFT;
#else
  FFT_storage FFT;
#endif
};


#define LOCALCE_MAXVDIM 2
struct union_cutoff { double constant, b, asqrtr, theor; };
struct union_cube { double constant, R, A, B, C, N, M, L; };
struct union_intrinsic { double A0, A2, B, MAX; };

struct localvariab {
  double R;
  int msg;
  int a;
  union {
    union_cutoff cutoff;
    union_cube cube;
    union_intrinsic intrinsic;
  };
};
typedef localvariab localvariabArray[LOCALCE_MAXVDIM * LOCALCE_MAXVDIM];
struct localCE_storage {
  double correction[MAXCEDIM * MAXCEDIM];
  localvariabArray q, q2;
};


struct approxCE_storage {
  int *idx;
};
Ulong NiceFFTNumber(Ulong nn);
int fastfourierInit(int *m, int dim, FFT_storage *FFT);
int fastfourier(double *data, int *m, int dim, bool inverse, FFT_storage *S);
int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *FFT);
void FFT_destruct(FFT_storage *S);
void FFT_NULL(FFT_storage *S);


struct trend_storage {
  double *x;
  int *xi;
  double *evalplane;
  int *powmatrix;
};

// see tbm.cc			  
struct tbm_storage {
  // bool genuine_dim[MAX362DIM];
  int ce_dim, simuspatialdim, method, spatialtotalpts, err;
  double center[MAXTBMSPDIM],  
    linesimuscale, linesimufactor, xline_length; 
};

// see spectral.cc

void metropolis(model *cov, gen_storage *S, double *x);
int search_metropolis(model *cov, gen_storage *S);
// see direct.cc


struct direct_storage {
  double *G;
};


// see sequential.cc
struct sequ_storage {
  int back, totpnts, spatialpnts, ntime, initial, delta_back_MuT;
  double *U11, *U22, *MuT,  *G,   *Cov21, *Inv22;
  double *res0;
};


// nugget
struct nugget_storage {
    //double sqrtnugget;
  bool spatialnugget, simuspatialnugget, simugrid;
  int total, *pos, *reduced_dim, *prod_dim, *index, *datapos;
  double *red_field;
};


// dummy version, hyperplane
struct hyper_storage{
  double rx[MAXHYPERDIM], center[MAXHYPERDIM], radius;
  hyper_pp_fct hyperplane;
};


struct plus_storage{
  bool conform[MAXSUB];
};


struct union_m3 {
  int **countvector, vertnumber, next_am_check;
  
  double minradius,  *lowerbounds, *areamatrix, *logvertnumber,
	*suppmin, *suppmax, radius,
    *loccentre; // only dummy variable! 	
};
struct union_shift {
  int  *mem2location, *loc2mem, *locindex, memcounter;
  double *location;
};  
struct union_normed {
  bool adaptive_nth, do_not_delete_C;
  int total, maxCi, nth, current_i, burnin, nCis;       
  Ulong zaehler, accepted;
  double *current_prob,  *current_cumprob, // *field,
    **C, *dummyCi, fmaxDfprop, max;
 };  
struct br_storage {
  int trendlen, zeropos;
  double **trend;
  int nr;

  union {    
    union_m3 m3;
    union_shift shift;
    union_normed normed;
  };

};


struct get_storage {
  int param_nr, size, vdim[2], // for nrow and ncol
    *idx;
  bool all;
};


struct pgs_storage {
  // urpsprunglich nur fuer pts_given_shape; jedoch allgemein
  // fuer shape functionen und zur Berechnung der Covariance/Variogram
  bool flathull, estimated_zhou_c, logmean; 
  double old_zhou, // for mcmc only
  totalmass, // (inverser) Normierungsfaktor, um Raeumliche Funktion
  //            zu einer Dichte zu machen
    currentthreshold, log_density, globalmin, intensity, alpha;
  int 
  rowscols, /// just for control
    size; //  *len, // global
 
  long double sum_zhou_c, sq_zhou_c;
  Long n_zhou_c;
  double zhou_c; // c in oesting, s, zhoy

  // Huetchen
  double *v, *y;  // local
  coord_type xgr;
  int *pos, // local
     *min, *max;// local dompp
  double *single, *total, *halfstepvector,  // global
    *localmin, *localmax, // local
    *minmean, *maxmean; // standard_shape
 
  // rf_interface.cc
   double
     *supportmin, *supportmax, *supportcentre, // global inkl. dompp
      *own_grid_start, *own_grid_step,
     *own_grid_len; // only for HUETCHENOWNGRIDSIZE
  int *gridlen, *end, *start, *delta, *nx, *own_grid_cumsum;
  double *xstart, *x, *inc;// local dompp
};



struct fctn_storage {
  // urpsprunglich nur fuer pts_given_shape; jedoch allgemein
  // fuer shape functionen und zur Berechnung der Covariance/Variogram
  int  rowscols,
    *cum,
    *end, *start, *nx,
    *endy, *startny, *ny;
 
  double *xstart, *x, *inc, *ystart, *y, *incy,
    *cross,  
    *C0x, *C0y;
  // param_set_fct param_set;

};


struct set_storage {
  param_set_fct set;
  //  void **valueRemote,
  //    **valueLocal;
  int variant;
  //     *bytes; // to be transfered
};

struct model_storage { // OK
  // Modelle nur hier, da diese rueckverfolgt werden muessen,
  // wenn prevloc umgesetzt wird
  model *cov, // OK, spezial, nie anders als  huetchen.cc verwenden
    *vario,
    *orig,
    *get_cov,
    *remote,
    *keys[MAXSUB];
  bool keys_given;
};

#define MAX_LIN_COMP (MAXSUB * MAXSUB)
#define model_undefined -1
#define model_morethan1 -2
typedef char NAname_type[MAX_NA][255];
struct likelihood_facts {
  int varmodel, NAs, nas[MAX_LIN_COMP],
    effect[MAX_LIN_COMP];
  model *Var; // ja nicht free, da nur pointer
  double *Matrix, 
    *pt_variance; // ja nicht free, da nur pointer
  bool trans_inv, isotropic, globalvariance;
  int newxdim, neffect;
  NAname_type names;
};

struct likelihood_storage {
  listoftype *datasets;  
  double **X, **YhatWithoutNA, *XCY, *XtX, *XitXi, *C, *CinvXY, *matrix,
    *betavec, **where_fixed, *sumY, *work, *Cwork,  *Xwork;
  int sets, fixedtrends, dettrends, random, max_total_data, *data_nas, maxbeta,
    cum_n_betas[MAX_LIN_COMP + 1], nas[MAX_LIN_COMP], nas_det[MAX_LIN_COMP], 
    nas_fixed[MAX_LIN_COMP], // # of additional NAs in a fixed effect
    nas_random[MAX_LIN_COMP],// # of additional NAs in a random effect
    nas_boxcox,
    nas_boxcox_mu; 
  bool dettrend_has_nas, fixedtrend_has_nas, random_has_nas, 
    data_has_nas,
    betas_separate, ignore_trend;
  char *betanames[MAX_LIN_COMP];
  model *fixed_effect[MAX_LIN_COMP], *det_effect[MAX_LIN_COMP], *rand_effect[MAX_LIN_COMP];
  likelihood_facts facts;
};



#define MAXDIM_POLY 2
struct polygon_storage {
  polygon *P;
  double **vdual;
  vertex *vprim;
  int n_vdual, n_vertex, n_v;
};


struct rect_storage {
  double inner, inner_const, inner_pow,
    outer, outer_const, outer_pow, outer_pow_const, step,
    *value, *weight, *tmp_weight, *right_endpoint, *ysort, *z;
  int nstep, tmp_n,
    *squeezed_dim, *asSign, *idx;
};


struct dollar_storage {
  bool busy, done, warned, timeprojection, separable,
    somegrid, blockmatrix;
  matrix_type type;
  double *sd, *save_aniso, *inv_aniso;
  int pid, *proj, nproj, 
    *cumsum, *total, *len, 
    n_z, n_z2, n_nx;
  isotropy_type orig_owniso;
  bool simplevar;
} ;

struct gatter_storage {
  double *z, *z1, *z2;
  int n_z, n_z1, n_z2;
};

struct earth_storage {
  double 
  P[9], cart_zenit[3]; // earth2cart u.ae.
};

struct extra_storage {
  double *a1, *a2, *a3,
    *b1, *b2,
    *c1, *c2;
  int *i1, //i2,
    *j1, //*j2,
    *k1, *k2,
    n_a1, n_a2, n_a3, n_b1, n_b2, n_c1, n_c2,
    n_i1, //n_i2,
    n_j1,// n_j2,
    n_k1, n_k2;
};

struct biwm_storage {
  bool nudiag_given, cdiag_given, donotcheck;// 13.11.2017 : usr_bool instead of bool
  double a[3], lg[3], aorig[3], nunew[3],
    scale[4], gamma[4], c[4];
};



struct bistable_storage {
  bool alphadiag_given, rhored_given;
  double alpha[3], scale[3], cdiag[2], rho, rhomax, rhored;
};


struct scatter_storage{
  int vdim, dim, *min, *max;
  double *step;
};


struct mcmc_storage{
  //bool done;
  //  int ;
  double integral, posvalue, *pos, *deltapos, *proposed, *propdelta;
};



struct mpp_properties {
  double 
    unnormedmass, // RRrectangle: mass of function that is at least
  //                as large as the given unnormed function
  // the equation
  //        unnormedmass * maxheight(of normed fctn) = mM[0]
  // frequently holds (in particular if mM[0]=1), but not always
  // exceptions is Loc where the unnormedmass is not calibrated against
  //        max(fct) = 1
  // calibration against the corresponding shape function does not help
  // (it would in case of Loc and Power$), but is does not in Rect,
  // since the later only guarantees that it is above the given function
  // (although ideally only slightly above or equal)
    
   // or of f / g ( Oesting, Sch    lather, Zhou), depending on the 
  // function (SHAPE_FCT)
    maxheights[MAXMPPVDIM], // maximum of f resp. of \d F // maxv=min(vdim, max)
    *mM, *mMplus // = int f^k \D \lambda, falls keine Verteilungsfamilie
    //                      und falls f kein stochastischer Prozess
    //       = \EE X^k = \int x^k \F(x) falls Verteilungsfamilie oder falls
    //             \EE X(0)^k, falls stochastischer Prozess
    // refradius , refsd,
    // totalmass 
    ;
  int // methnr,
    moments;
  //bool loc_done;
};


struct gen_storage { // cov_storage, do_storage, init_storage
  // wird in Struct initialisiert, sofern INIT aufgerufen wird;
  // abgespeichert wird es immer im aktuell aufrufenden (fuer next oder key)
  bool check, dosimulate; // used in biWM, BiGneiting
  spectral_storage Sspectral; // used in do
};



struct covariate_storage { // cov_storage, do_storage, init_storage
  // wird in Struct initialisiert, sofern INIT aufgerufen wird;
  // abgespeichert wird es immer im aktuell aufrufenden (fuer next oder key)
  location_type **loc;
  double *x;
  int matrix_err,
    pts; // used in var2cov
  int *nrow;  
  bool merged_x,
    nrow_checked,
    onlyOne,
    coordinate_mismatch;
};



struct bubble_storage {
  double *tau;
  int *rank, *start, *end;
};



///////////////////////////////////////////////////////////////////////
// AVL
#define AVL_FUNC_TYPES 1
struct cell_type {
    unsigned int *code;
    double colour;
};
typedef int (*avl_comparison_func) (cell_type *, cell_type *, int *);
typedef void (*avl_node_func) (cell_type*, int *);
typedef cell_type *(*avl_copy_func) (cell_type*, int *);


typedef
struct KEY_type KEY_type;
struct KEY_type {
  KEY_type *next;
  globalparam global;
  utilsparam global_utils;
  errorloc_type error_location;

  model *KEY[MODEL_MAX + 1];
  int pid, currentRegister, visitingpid, nzero, zaehler, set,
    use_external_boxcox;
  bool ok, stored_init,
    naok_range; // default =false;
  raw_type rawConcerns;
  char PREF_FAILURE[90 * Nothing];
  double *zerox;
  model *error_causing_cov;
};
extern KEY_type *PIDKEY[PIDMODULUS];
model **KEY();
KEY_type *KEYT();
int currentRegister();
void set_currentRegister(int cR);




// MLE.cc, gausdlikeli.cc
int internalSetAndGet(model *cov, int shortlen, 
		      int allowforintegerNA, bool excludetrend, // IN
		      int newxdim,  //IN
		      usr_bool globvar,
		      likelihood_facts *facts, // OUT	
		      sort_origin original);


#define booleanRange(IDX)\
  range->pmin[IDX] = range->min[IDX] = 0.0;	\
  range->pmax[IDX] = range->max[IDX] = 1.0;	\
  range->openmax[IDX] = range->openmin[IDX] = false


#define TH(i) (i) == 1 ? "st" : (i) == 2 ? "nd" : (i) == 3 ? "rd" : "th"

// own kappa name
#define OWNKAPPA(C, i)				\
  (STRCMP(C->kappanames[i], FREEVARIABLE) ? C->kappanames[i]		\
       : cov->ownkappanames != NULL && cov->ownkappanames[i] != NULL	\
       ? cov->ownkappanames[i] : "")



/* ********************************************************************** */
/*                           Functions                                    */
/* ********************************************************************** */


/* function necessary to set CircEmbed correctly if the function is not 
   even in some coordinates! 
*/
void InitModelList();   
/* absolutely necessary to call this function at the very beginning !
   but done automatically if SimulateRF, etc. is called for the first time
   (direct call necessary if further covariance functions ate to be added by 
   user!)
   initiating DefList with a dozen standard covariance models
*/


//void ExitInit(int err, bool *);
//void EnterInit(char *name);


#define GRIDEXPAND_AVOID Nan // cf TransformLoc
// to do variogram.cc:  statt True nur noch Avoid

#define DOLLAR_IMPOSSIBLE Nan

// Brown.cc, extremes.cc, nugget.cc, operator.cc, RMS.cc, tbm.cc, 
void TransformLoc(model *cov, bool timesep, usr_bool gridexpand,
		    usr_bool involvedollar);
// sequential.cc, variogram.cc
int TransformLoc(model *cov, double **xx, usr_bool involvedollar);
// variogram.h, operateor.cc 
int TransformLoc(model *cov, location_type *loc, double **xx, double **yy, usr_bool involvedollar);
// extremes
void TransformLocXY(model *cov, bool timesep, usr_bool gridexpand, 
		    usr_bool involvedollar);


listoftype *LIST_CREATE(int len, int type);
void LIST_DELETE(listoftype **x);
void listcpy(listoftype **To, listoftype *p, bool force_allocating);
void paramcpy(model *current, model *cov, bool freeing,
	      bool allocating, bool copy_lists, bool recursive, bool copy_mpp);
int covcpy(model **localcov, model *cov);
int covcpy(model **localcov, model *cov, bool copy_lists);
//int covcpy(model **localcov, bool sub, model *cov, location_type **prevloc);
int covcpy(model **localcov, bool sub, model *cov, // err
	   location_type **prevloc, location_type **ownloc,
	   bool copy_lists,  bool copy_randomparam, bool allowCopyingInterface);
int covcpyX(model **localcov, model *cov, // rename to covcpy!!
	   double *x, double *T, int spatialdim, int xdim, Long spatialpoints, bool Time, 
	   bool grid, bool distances);
int covcpyWithoutRandomParam(model **localcov, model *cov);
void Ssetcpy(model *localcov, model *remotecov, model *cov,
	     model *rmt);

void TaylorCopy(model *to, model *from);

int getmodelnr(char *name);
void setdefault(model *cov, int vdim0, int vdim1);
void setbackward(model *cov, model *sub);

int getmodelnr(char *name);

int setgrid(coord_type xgr, double *x, int spatialdim);
int partial_loc_set(location_type *loc, double *x, double *y,
		    Long spatialpoints, Long spatialpointsY, bool dist, int xdimOZ,
		    double *T, double *Ty,
		    bool grid, bool gridy, bool cpy);

// main fctn: used by userinferfaces.cc
int loc_set(double *x, double *y, double *T, double *Ty,
	    int spatialdim, // spatial dim only ! 
	    int xdimOZ,  Long spatialpoints, Long spatialpointsY,
	    bool Time, bool grid, bool gridY,
	    bool distances,
	    location_type **Loc);

// used by RMS
int loc_set(double *x, double *y, double *T, double *Ty,
	    int spatialdim, // spatial dim only ! 
	    int xdimOZ, Long spatialpoints, Long spatialpointsY,
	    bool Time, bool grid, bool gridY,
	    bool distances,
	    model *cov);

// used by Brown.cc, tbm.cc,
// getNset.cc : TransformCovLoc, 
int loc_set(double *x, double *T, 
	    int spatialdim, // spatial dim only !
	    int xdimOZ, // original ! 
	    Long spatialpoints, bool Time, bool grid,
	    bool distances,
	    model *cov);

// used by predict
void loc_set(SEXP ylist, location_type **Loc);

// used by operator.cc, shape.cc, userinterfaces.cc
location_type ** loc_set(SEXP xlist); 
location_type ** loc_set(location_type **Loc);
int empty_loc_set(model *cov, int dim, Long totalpoints, Long totalpointsY);

  
/*
int loc_set(double *x, double *T, 
	    int spatialdim, // spatial dim only !
	    int xdimOZ, // original ! 
	    Long spatialpoints, bool Time, bool grid,
	    bool distances, int n,
	    location_type ***Loc);
*/


void CheckModel(SEXP Model, double *x, double *y, double *T,  double *Ty,
		int spatialdim, /* spatial dim only ! */
		int xdim, int spatialpoints,  int spatialpointsY,
		bool grid, bool gridy,
		bool distances, bool Time, 
		SEXP xlist,
		raw_type raw,
		KEY_type *KT, int reg);

matrix_type Type(double *m, int nrow, int ncol);
double GetDiameter(location_type *loc);
double GetDiameter(location_type *loc, double *min, double *max,double *center);
double GetDiameter(location_type *loc, double *min, double *max,
		   double *center, bool docaniso, bool center_on_loc,
		   int *position);


void addModelX(model **pcov, int covnr);
void addModel(model *pcov, int subnr, int covnr);
void addModelKey(model *cov, int covnr);
void addModelKappa(model *pcov, int subnr, int covnr);
void addModel(model **pcov, int covnr, model *calling);
void addModel(model **pcov, int covnr, model *calling, bool nullOK); 
int addUnifModel(model *cov, double radius, model **newmodel);

void addVariable(char *name, double *x, int nrow, int ncol, SEXP env);
void addIntVariable(char *name, int *x, int nrow, int ncol, SEXP env);


void E1(spectral_storage *s, double A, double *e);
void E2(spectral_storage *s, double A, double *e);
void E12(spectral_storage *s, int dim, double A, double *e);
void E3(spectral_storage *s, double A, double *e);
void E(int dim, spectral_storage *s, double A, double *e);



void addmsg(double value, const char *Sign, double y, char *msg);
int checkkappas(model *cov);
int checkkappas(model *cov, bool on);

//double gaussInt(int d, int xi, double sigma, double R);

void updatepref(model *cov, model *sub);

void leer(int level);

int get_ranges(model *cov, model **min, model **max,
		model **pmin, model **pmax, 
		model **openmin, model **openmax);
int check_recursive_range(model *cov, bool NAOK);

void xtime2x(double *x, int nx, double *T, int len, double **newx, int nrow);


//int Checking(model **Cov);
int check2passTF(model *cov, system_type* s, Types type, int vdim, Types frame);
int check2passframe(model *cov, system_type* s, int vdim0, int vdim1, 
		 Types frame);
int check2passtype(model *cov, system_type* s, Types type,
		   int vdim0, int vdim1, Types frame);
int check2Xnotrafo(model *cov, int logicaldim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
		   int vdim, Types frame);
int check2X(model *cov, int logicaldim, int tsxdim, Types type, 
	    domain_type domprev, isotropy_type isoprev, 
	    int vdim, Types frame);
int check2X(model *cov, int logicaldim, int tsxdim,
	    Types type, domain_type domprev, isotropy_type isoprev,
	    int *vdim, Types frame);
int check2X(model *cov, int logicaldim, int tsxdim, Types type, 
	    domain_type domprev, isotropy_type isoprev,
	    int vdim0, int vdim1, Types frame, bool coordinate_trafo);
int check2X(model *cov, int vdim0, int vdim1, Types frame, bool coord_trafo);
int check2Xthroughout(model *cov, model *prev, 
	   Types type, domain_type domprev, isotropy_type isoprev,
		  int vdim, Types frame);
int CheckPosVar(model *cov, Types type, domain_type dom, isotropy_type iso,
		int vdim, Types frame);


#define CHECK_GEN(C,V0, V1, F, CT) check2X(C,V0, V1, F, CT)
#define CHECK_ONLY(C) CHECK_GEN(C,(C)->vdim[0],(C)->vdim[1],(C)->frame,false)
#define CHECK_NOPASS(C) check2passframe(C, OWN, VDIM0, VDIM1, cov->frame)
#define CHECK_PASSTF(C,T,V,F) check2passTF(C, OWN, T, V, F)
#define CHECK_PASSFRAME(C,F) check2passframe(C, OWN, VDIM0, VDIM1, F)
#define CHECK_PASSTYPE(C,T) check2passtype(C, OWN, T, VDIM0, VDIM1, cov->frame)
#define CHECK_VDIM(C,T,X,type,D,I,V0,V1,F) check2X(C,T,X,type,D,I,V0,V1,F,true)
#define CHECKPOSVAR(C,type,D,I,V,F) CheckPosVar(C,type,D,I,V,F)
#define CHECK_R(C, vdim)				     \
  CHECK_VDIM(C, vdim, vdim, RandomType, KERNEL, CARTESIAN_COORD, \
	  vdim, 1, RandomType)			       


int INIT_RANDOM_intern(model *M, int moments, gen_storage *s, double *p);
int INIT_intern(model *M, int moments, gen_storage *s);
int REINIT_intern(model *M, int moments, gen_storage *s);
 

void kdefault(model *cov, int i, double v);

// KeyFacts
void iexplDollar(model *cov, bool MLEnatsc_only);


// allocs
int alloc_mpp_M(model *cov, int moments);
void free_mpp_M(model *cov);
int alloc_pgs(model *cov);
int alloc_fctn(model *cov, int dim, int rowscols);

// others
void crash();


double *getAnisoMatrix(model *cov, int *nrow, int *ncol);
double *getAnisoMatrix(model *cov, bool null_if_id, int *nrow, int *ncol);


void SetLoc2NewLoc(model *cov, location_type **Loc);
void SetLoc2NewLoc(model *cov, location_type **old, location_type **neu,
		   int alle);

int ReturnOwnField(model *cov);
int ReturnOtherField(model *cov, model *which);



///////////////////////////////////////////////////////////////////////
// UNAUFGERAEUMT:
///////////////////////////////////////////////////////////////////////




#define MAXDEFMATRIX 3



#define AveMaxDim 10 
#define CoxMaxDim 3 
#define StpMaxDim 10 
#define EaxxaMaxDim 10  
#define ShiftMaxDim 10  
#define ParsWMMaxVDim 10 

  
struct model {
  // user given information
  int err_level, err, 
    zaehler; /* for debugging only */
  errorstring_type err_msg;  
  param_type px;       // 24 b
  int nrow[MAXPARAM], // 24 bytes
    ncol[MAXPARAM];   // 24 bytes
  double *q;
  int qlen, variant,
    nsub; /* current number of submodels */
  model *sub[MAXSUB], *kappasub[MAXPARAM], *calling, *root, *key; 
  KEY_type *base;
  location_type **prevloc, **ownloc;

  char **ownkappanames;
  user_given_type user_given;
  Systems_type prev, /* by previous model: logdim: required
			maxdim: NA
			xdim  : intended to be delivered
			cumxdim: (calculated)
			type  : required
			domain: to be delivered/required
			iso   : to be delivered/required 
			nr: refers to trafos among coordinate systems, if any
			    else UNSET

		     */
    gatter, /* copy of prev if no coordinate transfomation necessary
             otherwise the new coordinate system is given
             note that for gatter, only last, xdim, and logicaldim
	     are defined. Other values are arbitrary
	     nr: refers to trafos within coordinate systems; 
	    */
    own; /* 
	    status after all gatter-calls, as required by one of the
	    definitions of the model 
	    nr : contains ones own covariance nr
	 */

  /////////////////////////////////////////////////////////
  // VARIABLES PASSED DOWNWARDS
  /////////////////////////////////////////////////////////
  Types frame; /* current frame of reference/surrounding environment
		  of the model, determined by the calling model */

  int
     vdim[2],
     full_derivs, rese_derivs; /* rechtsseitige bzw. vollstaendige Anzahl 
				 Ableitungen number of derivatives of model 
				 combined with information of the submodels */

  monotone_type
     monotone;/* for simple model: normal mix model iff maxdim = INFDIM
		     for hypermodel we need nonetheless this parameter
		     parameter set by getinfo(), not within info of the cov 
		     fcts comining information of submodels  
	       */
 

  /////////////////////////////////////////////////////////
  // VARIABLES PASSED UPWARDS
  /////////////////////////////////////////////////////////
  // forward analysis of user's information
 
  double  *rf, // for storing random shapes; this is a pointer pointing
  //               to *rf in storage (logically rf only belongs to 
  //               a single atom. But programmed as such, this would take
  //               an enormous amount of memory in case of M3 processes)
    logspeed; /* 
		      logspeed = lim_{h->infty} \gamma(h)/l o g(h) 
		      in case of isotropic model and RF_NAN otherwise
		   */

  ext_bool fieldreturn,
    finiterange, /* also information obtained by model and submodels 
		  */
    loggiven;
  
  bool
    randomkappa, // 
    matrix_indep_of_x,
    hess,  /* can a hessian matrix be provided? */
    initialised, // is the simulation initialised? 
    origrf,      // does *rf point to allocated memory?
    checked;     // has the model been checked successfully?

  allowedD_type allowedD;
  allowedI_type allowedI;
  bool IallowedDone, DallowedDone;
   
  double taylor[MAXTAYLOR][TaylorPow + 1], 
    tail[MAXTAYLOR][TaylorExpPow + 1]; 
  int taylorN, // number of summands in the taylor expansion -- whatever the
  // exponents of the terms are
    tailN; // dito


  pref_type pref; /* 0 : not possible; 
		     5 : best possible
		     (including subs)!
		   */
  //   user; /* 0 : not possible; 
  //	     5 : best possible
  //	     (including subs)!
  //		  */

  // Information used and created within simuations
  Methods method; /* the current method (out of SimulationsMeth) which 
			    is tried at the moment or which has been 
			    successfully initialized */
  mpp_properties mpp; 
  ptwise_type ptwise_definite;

  simu_storage simu;
  ce_storage *Sce;
  localCE_storage *SlocalCE;
  approxCE_storage *SapproxCE;
  direct_storage *Sdirect;
  hyper_storage *Shyper;
  nugget_storage *Snugget;
  plus_storage *Splus;
  sequ_storage *Ssequ;
  tbm_storage *Stbm;
   br_storage *Sbr;
  get_storage *Sget;
  pgs_storage *Spgs;
  fctn_storage *Sfctn;
  set_storage *Sset;
  polygon_storage *Spolygon;
  rect_storage *Srect;
  dollar_storage *Sdollar;
  gatter_storage *Sgatter;
  earth_storage *Searth;
  extra_storage *Sextra;
  solve_storage *Ssolve;
  biwm_storage *Sbiwm;
  bistable_storage *Sbistable;
  scatter_storage *Sscatter;
  mcmc_storage *Smcmc;
  likelihood_storage *Slikelihood;
  covariate_storage *Scovariate; // only once for the whole model tree
  bubble_storage *Sbubble;
  mle_storage *Smle;
  model_storage *Smodel;
  trend_storage *Strend;
  gen_storage *Sgen; // only once for the whole model tree

  //check in getNset:  COV_DELETE_WITHOUTSUB, COV_ALWAYS_NULL; add *_DEL, *_NULL
  
  //select_storage *Sselect;
};


void KEY_type_NULL(KEY_type  *x);
void KEY_type_DELETE(KEY_type **S);
void LOC_DELETE(location_type ***Loc); // OK
void LOC_SINGLE_NULL(location_type *loc, int len, int dim);
//location_type **LOCLIST_CREATE(int n, int dim);
void LOC_SINGLE_DELETE(location_type **Loc);
void LOC_SINGLE_DELETE_Y(location_type *loc);
 location_type **LOCLIST_CREATE(int n, int dim);
void COV_ALWAYS_NULL(model *cov);
void SYSTEM_NULL(system_type *sys, int len);
void COV_DELETE_(model **cov, model *save);
void COV_NULL(model *cov, KEY_type *base);
void COV_DELETE_WITHOUTSUB(model **Cov, model *save);
void COV_DELETE_WITHOUT_LOC(model **Cov, model *save);
void ce_NULL(ce_storage* x);
void ce_DELETE(ce_storage **S);
void localCE_NULL(localCE_storage* x);
void localCE_DELETE(localCE_storage**S);
void approxCE_NULL(approxCE_storage* x);
void approxCE_DELETE(approxCE_storage **S);
void direct_NULL(direct_storage  *x);
void direct_DELETE(direct_storage  ** S);
void hyper_NULL(hyper_storage* x);
void hyper_DELETE(hyper_storage  **S);
void nugget_NULL(nugget_storage *x);
void nugget_DELETE(nugget_storage ** S);
void plus_NULL(plus_storage *x);
void plus_DELETE(plus_storage ** S);
void sequ_NULL(sequ_storage *x);
void sequ_DELETE(sequ_storage **S);
void spectral_NULL(spectral_storage *x);
void spectral_DELETE(spectral_storage **S);
void tbm_DELETE(tbm_storage **S); 
void tbm_NULL(tbm_storage* x);
void br_DELETE(br_storage **S); 
void br_NULL(br_storage* x);
void get_NULL(get_storage *S);
void get_DELETE(get_storage **S);
void pgs_DELETE(pgs_storage **S); 
void pgs_NULL(pgs_storage* x);
void fctn_DELETE(fctn_storage **S); 
void fctn_NULL(fctn_storage* x);
void set_DELETE(set_storage **S); 
void set_NULL(set_storage* x);
void polygon_DELETE(polygon_storage **S); 
void polygon_NULL(polygon_storage* x);
void rect_DELETE(rect_storage **S); 
void rect_NULL(rect_storage* x);
void dollar_DELETE(dollar_storage **S); 
void dollar_NULL(dollar_storage* x);
void gatter_DELETE(gatter_storage **S); 
void gatter_NULL(gatter_storage* x);
void earth_DELETE(earth_storage **S); 
void earth_NULL(earth_storage* x);
void extra_DELETE(extra_storage **S); 
void extra_NULL(extra_storage* x);
void biwm_DELETE(biwm_storage **S); 
void biwm_NULL(biwm_storage* x);
void bistable_DELETE(bistable_storage **S);
void bistable_NULL(bistable_storage* x);
void scatter_DELETE(scatter_storage **S);
void scatter_NULL(scatter_storage* x);
void mcmc_DELETE(mcmc_storage **S);
void mcmc_NULL(mcmc_storage* x);
void likelihood_NULL(likelihood_storage *x);
void likelihood_DELETE(likelihood_storage **S);
//likelihood_storage * likelihood_CREATE(int n);
void likelihood_facts_NULL(likelihood_facts *x);
void likelihood_facts_DELETE(likelihood_facts *x);
void covariate_NULL(covariate_storage *x);
void covariate_DELETE(covariate_storage **S);
void bubble_NULL(bubble_storage *x);
void bubble_DELETE(bubble_storage **S);
void mle_NULL(mle_storage *x);
void mle_DELETE(mle_storage **S);
void model_NULL(model_storage *x);
void model_DELETE(model_storage **S, model *save);
void trend_DELETE(trend_storage ** S);
void trend_NULL(trend_storage* x);
void gen_NULL(gen_storage *x);
void gen_DELETE(gen_storage **S);
// grep Strend *; grep trend_NULL *

void BRtrend_destruct(double **BRtrend, int trendlen);
int StructBR(model *cov, gen_storage *s, model **atom);
void FFT_destruct(FFT_storage *FFT);


double *ZERO(model *cov);
double *ZERO(int dim, KEY_type *KT);

//void SELECT_NULL(select_storage *S);
//void SELECT_DELETE(select_storage **S);


#define Nick(Cov) (DefList[MODELNR(Cov)].nick)
//#define NICK(Cov) (isDollar(Cov) ? Nick((Cov)->sub[0]) : isPlusMal(Cov)  ? NAME(Cov) : Nick(Cov))
#define NICK(Cov) (isDollar(Cov) ? Nick((Cov)->sub[0]) : Nick(Cov))
#define NAME(Cov) DefList[MODELNR(Cov)].name
#define KNAME(NAME) DefList[COVNR].kappanames[NAME]
#define SNAME(NAME) DefList[COVNR].subnames[NAME]

int rPoissonPolygon2(polygon_storage *S, double lambda, bool do_centering);

//////////////////////////////////////////////////////////////////////
// DEDUGGING INFORMATION
//////////////////////////////////////////////////////////////////////


///////////////////////////////////////////
#define SET_DESTRUCT(A)\
  assert(meth->S==NULL && meth->destruct==NULL);	\
  meth->destruct = A;
								

///////////////////////////////////////////


#define COV_DELETE COV_DELETE_
#define INIT(Cov, Moments, S) INIT_intern(Cov, Moments, S)
#define REINIT(Cov, Moments, S) REINIT_intern(Cov, Moments, S)
#define INIT_RANDOM(Cov, Moments, S, P) INIT_RANDOM_intern(Cov, Moments, S, P)




#define DO(Cov, S) {						\
    assert((Cov)->initialised);				\
    ASSERT_GATTER(Cov);							\
    PL--;								\
    DefList[FIRSTGATTER].Do(Cov, S);					\
    PL++;								\
  }
   
#define DORANDOM(Cov, S) {						\
    assert(!isBad(TypeConsistency(RandomType, Cov, ISO(SYSOF(Cov), 0)))); \
    ASSERT_GATTER(Cov);							\
    PL--;								\
    DefList[FIRSTGATTER].DoRandom(Cov, S);				\
    PL++;								\
  }

#define COV(X, Info, Cov, V) {						\
    ASSERT_GATTER(Cov); ASSERT_XONLY(Cov); DefList[FIRSTGATTER].cov(X, Info, Cov, V);}
#define LOGCOV(X, Info, Cov, V, S) {				\
    ASSERT_GATTER(Cov);ASSERT_XONLY(Cov);DefList[FIRSTGATTER].log(X, Info, Cov, V, S);}
#define SHAPE COV
#define FCTN COV
#define ABSFCTN(X, Info, Cov, V) {ASSERT_XONLY(Cov); COV(X, Info, Cov, V); *(V) = FABS(*(V)); }
#define LOGSHAPE LOGCOV

void nonstat2statcov(double *x, double *y, int*info, model *cov, double *v);
void lognonstat2statcov(double *x, double *y, int*info, model *cov, double *v,
		   double *sign);
#define NONSTAT2STATCOV(X, Y, Info, Cov, V) {ASSERT_XONLY(Cov);nonstat2statcov(X, Y, Info, Cov, V);}
#define LOGNONSTAT2STATCOV(X, Y, Info, Cov, V, S) {ASSERT_XONLY(Cov);	\
    lognonstat2statcov(X, Y, Info, Cov, V, S);}

#define NONSTATCOV(X, Y, Info, Cov, V) {ASSERT_KERNEL(Cov);		\
    ASSERT_GATTER(Cov); DefList[FIRSTGATTER].nonstat_cov(X, Y, Info, Cov,V);}
#define LOGNONSTATCOV(X, Y, Info, Cov, V, S) {	ASSERT_KERNEL(Cov);	\
    ASSERT_GATTER(Cov);DefList[FIRSTGATTER].nonstatlog(X, Y, Info, Cov,V,S);}
#define Abl1(X, Info, Cov, V) {ASSERT_XONLY(Cov); \
    ASSERT_GATTER(Cov);DefList[FIRSTGATTER].D(X, Info, Cov, V);}
#define Abl2(X, Info, Cov, V) {ASSERT_XONLY(Cov); \
    ASSERT_GATTER(Cov);DefList[FIRSTGATTER].D2(X, Info, Cov, V);}
#define Abl3(X, Info, Cov, V) {ASSERT_XONLY(Cov); \
    ASSERT_GATTER(Cov);DefList[FIRSTGATTER].D3(X, Info, Cov, V);} /* OK ? */
#define Abl4(X, Info, Cov, V) {ASSERT_XONLY(Cov); \
    ASSERT_GATTER(Cov);DefList[FIRSTGATTER].D4(X, Info, Cov, V);} /* OK ? */

#define INVERSE(V, Cov, X) {						\
    ASSERT_GATTER(Cov);DefList[FIRSTGATTER].inverse(V, Cov, X);}
#define INVERSENONSTAT(V, Cov, X, Y) {\
    ASSERT_GATTER(Cov); DefList[FIRSTGATTER].inverse_nonstat(V, Cov, X, Y);}
#define INVERSENONSTATLOG(V, Cov, X, Y) {\
    ASSERT_GATTER(Cov); DefList[FIRSTGATTER].nonstat_loginverse(V, Cov, X, Y);}

#define VTLG_D(X, Info, Cov, V) { /* ASSERT_XONLY(Cov); */		\
    ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].D(X, Info, Cov, V);}//kein gatter notw.
#define VTLG_DLOG(X, Info, Cov, V) { /* TO DO bei allen Verteilungsaufruefn statt ASSERT_XONLY(Cov); -- irgendeeine sicherheit waere gut*/		\
    ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].logD(X, Info, Cov, V);}
#define VTLG_P(X, Info, Cov, V) { /* ASSERT_XONLY(Cov); */		\
    ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].cov(X, Info, Cov, V);} 
#define VTLG_P2SIDED(X, Y, Info, Cov, V) {/* nicht gatter, da X=NULL sein kann !*/  \
    ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].nonstat_cov(X, Y, Info, Cov, V);} 
#define VTLG_Q(V, Cov, X) { \
    ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].inverse(V, Cov, X);}
#define VTLG_R(X, Info, Cov, V) { /* ASSERT_XONLY(Cov); */		\
    ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].random(X, Info, Cov, V);} /* dito */
#define VTLG_R2SIDED(X, Y, Info, Cov, V) {ASSERT_KERNEL(Cov);\
    ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].nonstat_random(X, Y, Info, Cov, V);}
#define INVERSENONSTAT_D(V, Cov, X, Y)		       \
  { ASSERT_CHECKED(Cov); DefList[MODELNR(Cov)].inverse_nonstat_D(V, Cov, X, Y);}


//#define DENSITYFCT(X, Dens, V) DefList[FIRSTGATTER].approx_dens(X, Dens, V)
#define SPECTRAL(Cov, S, E) {\
    ASSERT_GATTER(Cov); DefList[MODELNR(Cov)].spectral(Cov, S,E);}/*not gatter*/
#define TBM2CALL(X, Info, Cov, V) {ASSERT_KERNEL(Cov);			\
    ASSERT_GATTER(Cov); assert(DefList[MODELNR(Cov)].tbm2 != NULL);	\
    DefList[MODELNR(Cov)].tbm2(X, Info, Cov, V);}
#define HESSE(X, Info, Cov, V) {ASSERT_XONLY(Cov); \
    ASSERT_GATTER(Cov);DefList[MODELNR(Cov)].hess(X, Info, Cov, V);}
#define NABLA(X, Info, Cov, V) {ASSERT_XONLY(Cov); \
    ASSERT_GATTER(Cov);DefList[MODELNR(Cov)].nabla(X, Info, Cov, V);}

#define ASSERT_XONLY(Cov) assert(isXonly(PREVSYSOF(Cov)))
#define ASSERT_KERNEL(Cov) assert(isKernel(PREVSYSOF(Cov)))


//#define FRAME_ASSERT(F)				
//  if (has##F##Frame(cov)){ /* NICHT! : '(Frame)'*/	
//  } else {								
//    assert(({PMI(cov) ; true;}));					
//    SERR2("Frame '%.50s' not recognised by '%.50s'.",			
//	  TYPE_NAMES[cov->frame], NICK(cov));				
//  }

#define ILLEGAL_FRAME							\
  NERR4(ERRORILLEGALFRAME, "cannot initiate '%.50s' within frame '%.50s' [debug info: '%.50s' at line %d]", \
	NICK(cov), TYPE_NAMES[cov->frame], __FILE__, __LINE__)

#define ILLEGAL_FRAME_STRUCT \
  NERR2(ERRORILLEGALFRAME, "cannot restructure '%.50s' by frame '%.50s'",\
	NICK(cov), TYPE_NAMES[cov->frame])

#define ASSERT_NEWMODEL_NOT_NULL\
  if (newmodel != NULL) { } else 				\
    SERR1("unexpected call of struct_%.50s", NAME(cov));	       
  

#define ASSERT_NEWMODEL_NULL\
  if (newmodel == NULL) { } else 			\
    SERR1("Unexpected call of struct_%.50s", NAME(cov));
 

#define ASSERT_ONE_SUBMODEL(Cov) {					\
    { DEBUGINFO; }							\
    if (!(((Cov)->sub[0] == NULL) xor ((Cov)->sub[1] == NULL))) {	\
      defn *C = DefList + MODELNR(Cov);					\
      SERR2("either '%.50s' or '%.50s' must be given", C->subnames[0],	\
	    C->subnames[1]);						\
    }									\
  }

#define ASSERT_FRAME_DEFINED(Cov) { \
    DEBUGINFO;							\
    if (isBad((Cov)->frame))					\
      SERR1("'%.50s' has badly defined frame.", NICK(Cov));	\
  }


#define EXT_NEW_COV_STORAGE(cov, new) {				\
    assert((cov)->S##new == NULL);				\
   if ((cov)->S##new != NULL) {				\
      Ext_##new##_DELETE(&((cov)->S##new));		    	\
      assert((cov)->S##new == NULL);					\
    }									\
    (cov)->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
      Ext_##new##_NULL((cov)->S##new);					\
    if ((cov)->S##new == NULL) BUG;					\
  }					


#define EXT_NEW_STORAGE(new)	\
  EXT_NEW_COV_STORAGE(cov, new)

#define FORCE_NEW_COV_STORAGE(cov, new) {				\
    if ((cov)->S##new != NULL) { /* for safety */			\
      new##_DELETE(&((cov)->S##new));					\
      assert((cov)->S##new == NULL);					\
    }									\
    (cov)->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    new##_NULL((cov)->S##new);						\
    if ((cov)->S##new == NULL) BUG;					\
  }   	
#define FORCE_NEW_STORAGE(new)	\
  FORCE_NEW_COV_STORAGE(cov, new)


#define NEW_COV_STORAGE(cov, new) {					\
    assert((cov)->S##new == NULL);					\
    if ((cov)->S##new != NULL) { /* for safety */			\
      new##_DELETE(&((cov)->S##new));					\
      assert((cov)->S##new == NULL);					\
    }									\
    (cov)->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    new##_NULL((cov)->S##new);						\
    if ((cov)->S##new == NULL) BUG;					\
  }   	
#define NEW_STORAGE(new) NEW_COV_STORAGE(cov, new)


#define NEW_COV_STORAGE_SAVE(cov, new) {				\
    assert((cov)->S##new == NULL);					\
    if ((cov)->S##new != NULL) { /* should never be valid! */		\
      new##_DELETE(&((cov)->S##new), cov); /* unterschied !*/	\
      assert((cov)->S##new == NULL);					\
    }									\
    (cov)->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    new##_NULL((cov)->S##new);						\
    if ((cov)->S##new == NULL) BUG;					\
  }								


#define ONCE_NEW_COV_STORAGE(cov, new)					\
  if ((cov)->S##new == NULL) {						\
    (cov)->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    new##_NULL((cov)->S##new);						\
    if ((cov)->S##new == NULL) BUG;					\
  }								

#define ONCE_NEW_STORAGE(new)	\
  ONCE_NEW_COV_STORAGE(cov, new)


#define CONDCOV_NEW_STORAGE(cov, new, WHAT) {				\
  if (cov->S##new != NULL && cov->S##new->WHAT != NULL) {		\
    new##_DELETE(&(cov->S##new));					\
    assert((cov)->S##new == NULL);					\
  }								\
  if (cov->S##new == NULL) {					\
   cov->S##new = (new##_storage *) MALLOC(sizeof(new##_storage));	\
    new##_NULL(cov->S##new);					\
    if (cov->S##new == NULL) BUG;					\
  }								\
  assert((cov)->S##new->WHAT == NULL);				\
  }

#define COND_NEW_STORAGE(new, WHAT) CONDCOV_NEW_STORAGE(cov, new, WHAT)
#define SOLVE_STORAGE EXT_NEW_STORAGE(solve)

#define ALLCCOV_NEW(cov, Snew, Z, SIZE, WHAT)				\
  assert((cov)->Snew != NULL);						\
  double *Z = (cov)->Snew->WHAT;					\
  if ((Z) != NULL) { } else						\
    (Z) = (cov)->Snew->WHAT = (double*) MALLOC(sizeof(double) * (SIZE)) 
  
#define ALLC_NEW(Snew, Z, SIZE, WHAT) ALLCCOV_NEW(cov, Snew, Z, SIZE, WHAT)
  
#define ALLC_NEWINT(Snew, Z, SIZE, WHAT)				\
  assert((cov)->Snew != NULL);						\
   int *Z = cov->Snew->WHAT;						\
  if (Z != NULL) { } else						\
    Z = cov->Snew->WHAT = (int*) MALLOC(sizeof(int) * (SIZE))


#define DISTMAXSTEPS 1000
#define XSIZE 16
#define XXSIZE 116
#define XXXSIZE 1116 // >= DISTAMAXSTEPS + 2 *dim
#if (XXSIZE < 1L + MATERN_NU_THRES)
BUG
#endif



#define EXTRA_STORAGE NEW_STORAGE(extra)
#define ONCE_EXTRA_STORAGE ONCE_NEW_STORAGE(extra)
#define GATTER_STORAGE COND_NEW_STORAGE(gatter, z);
#define GATTER_STORAGE COND_NEW_STORAGE(gatter, z);
#define TALLOC_ASSERT(cov, Snew, SIZE, STANDARD)			\
   assert((cov)->Snew != NULL);						\
   assert((SIZE) > 0);							\
   assert((STANDARD) > 0)


#ifdef DO_TALLOC_SAVE
#define TALLOCCOV_NEW(cov,Snew, Z, SIZE, WHAT, STANDARD)	\
  TALLOC_ASSERT(cov, Snew, SIZE, STANDARD);				\
  bool free_##Z = ((SIZE) > (STANDARD));				\
  double WHAT##__X[STANDARD], *Z,  *WHAT##__Y = NULL;			\
  if (free_##Z) {							\
    WHAT##__Y = (double*) MALLOC(sizeof(double)*(SIZE));		\
    Z =  WHAT##__Y;							\
  } else Z = WHAT##__X

#define TALLOCCOV_G_NEW(cov, Z, SIZE, STANDARD)			\
  assert((SIZE) > 0);							\
  assert((STANDARD) > 0);						\
  free_##Z = ((SIZE) > (STANDARD));					\
  if (free_##Z) {							\
    Z##__Y=(double*) MALLOC(sizeof(double)*(SIZE));			\
    Z = Z##__Y;								\
  } else Z = Z##__X


#define TALLOC_NEWINT(Snew, Z, SIZE, WHAT, STANDARD)		\
  TALLOC_ASSERT(cov, Snew, SIZE, STANDARD);				\
  bool free_##Z = ((SIZE) > (STANDARD));				\
  int WHAT##__X[STANDARD], *WHAT##__Y = NULL, *Z;			\
  if (free_##Z) {							\
    WHAT##__Y =(int*) MALLOC(sizeof(int)*(SIZE));			\
    Z =  WHAT##__Y;							\
  } else Z = WHAT##__X

  
#define END_TALLOC_NEW(WHAT) FREE(WHAT##__Y)//very crucial that Z isn't freed!
//                                           see stat2(_Intern) !
#define FREE_TALLOC(Z) if (free_##Z) FREE(Z)

#define TALLOC_DOUBLE(Z) bool free_##Z=false; double *Z=NULL, Z##__X[XSIZE], *Z##__Y=NULL;
#define TALLOC_INT(Z) bool free_##Z=false; int *Z=NULL, Z##__X[XSIZE], *Z##__Y=NULL
// MALLOC ?
// MALLOC :




#else // NOT DO_TALLOC_SAVE
#define TALLOCCOV_NEW(cov, Snew, Z, SIZE, WHAT, STANDARD)	\
  TALLOC_ASSERT(cov, Snew, SIZE, STANDARD);				\
   double *Z = (cov)->Snew->WHAT;					\
  /* if ((cov)->Snew->n_##WHAT < SIZE) { FREE(Z); (cov)->Snew->n_##WHAT = SIZE; } */\
  if ((Z) != NULL) { } else						\
    (Z) = (cov)->Snew->WHAT = (double*) MALLOC(sizeof(double) * (SIZE)) 

#define TALLOCCOV_G_NEW(cov, Snew, Z, SIZE, WHAT, STANDARD)		\
  TALLOC_ASSERT(cov, Snew, SIZE, STANDARD);				\
  Z = (cov)->Snew->WHAT;					\
  /* if ((cov)->Snew->n_##WHAT < SIZE) { FREE(Z); (cov)->Snew->n_##WHAT = SIZE; } */\
  if ((Z) != NULL) { } else						\
    (Z) = (cov)->Snew->WHAT = (double*) MALLOC(sizeof(double) * (SIZE)) 
  
#define TALLOC_NEWINT(Snew, Z, SIZE, WHAT, STANDARD)		\
  TALLOC_ASSERT(cov, Snew, SIZE, STANDARD);				\
  int *Z = cov->Snew->WHAT;						\
  if (Z != NULL) { } else						\
    Z = cov->Snew->WHAT = (int*) MALLOC(sizeof(int) * (SIZE)) 

#define END_TALLOC_NEW(WHAT) 
#define FREE_TALLOC(Z) 
#define TALLOC_DOUBLE(z) double *z=NULL
#define TALLOC_INT(z) int *z=NULL
#endif // (NOT) DO_TALLOC_SAVE


#define TALLOC_NEW(Snew, Z, SIZE, WHAT, STANDARD)		\
  TALLOCCOV_NEW(cov, Snew, Z, SIZE, WHAT, STANDARD)
// ACHTUNG!!! NIE TALLOC_GLOBAL_XX1 OHNE TALLOC_DOUBLE ABZUANENDERN !!
#ifdef DO_TALLOC_SAVE
#define TALLOC_GLOBAL_X1(Z, SIZE) TALLOCCOV_G_NEW(cov, Z, SIZE, XSIZE) 
#define TALLOC_GLOBAL_X2(Z, SIZE) TALLOCCOV_G_NEW(cov, Z, SIZE, XSIZE) 
#define TALLOC_GLOBAL_X3(Z, SIZE) TALLOCCOV_G_NEW(cov, Z, SIZE, XSIZE)
#define TALLOC_GATTER_GLOBAL(Z,SIZE) TALLOCCOV_G_NEW(cov, Z, SIZE, XSIZE)
#else
#define TALLOC_GATTER_GLOBAL(Z,SIZE) TALLOCCOV_G_NEW(cov,Sgatter,Z,SIZE,Z,XSIZE)
// !!! ACHTUNG ! bei Verwendung von TALLOC auf Sextra basierend, darf
// bis zu END_TALLOC kein cov stehen, oder es muss sichergestellt werden
// dass die Fkt mit argument 'cov' nicht auch noch TALLOC aufruft!!
#define TALLOC_GLOBAL_X1(Z, SIZE) TALLOCCOV_G_NEW(cov, Sextra, Z,SIZE,a1,XSIZE) 
#define TALLOC_GLOBAL_X2(Z, SIZE) TALLOCCOV_G_NEW(cov, Sextra, Z,SIZE,a2,XSIZE) 
#define TALLOC_GLOBAL_X3(Z, SIZE) TALLOCCOV_G_NEW(cov, Sextra, Z, SIZE,a3,XSIZE)
#endif
// !!! ACHTUNG ! bei Verwendung von TALLOC auf Sextra basierend, darf
// bis zu END_TALLOC kein cov stehen, oder es muss sichergestellt werden
// dass die Fkt mit argument 'cov' nicht auch noch TALLOC aufruft!!
#define TALLOC_X1(Z, SIZE) TALLOC_NEW(Sextra, Z, SIZE, a1, XSIZE) 
#define TALLOC_X2(Z, SIZE) TALLOC_NEW(Sextra, Z, SIZE, a2, XSIZE) 
#define TALLOC_X3(Z, SIZE) TALLOC_NEW(Sextra, Z, SIZE, a3, XSIZE) 
#define TALLOC_XX1(Z, SIZE) TALLOC_NEW(Sextra, Z, SIZE, b1, XXSIZE)
#define TALLOC_XX2(Z, SIZE) TALLOC_NEW(Sextra, Z, SIZE, b2, XXSIZE)
#define TALLOC_XXX1(Z, SIZE) TALLOC_NEW(Sextra, Z, SIZE,c1,XXXSIZE)
#define TALLOC_XXX2(Z, SIZE) TALLOC_NEW(Sextra, Z, SIZE,c2,XXXSIZE)
#define TALLOC_L1(Z, SIZE) TALLOC_NEWINT(Sextra, Z, SIZE, i1, XSIZE)
#define TALLOC_L2(Z, SIZE) TALLOC_NEWINT(Sextra, Z, SIZE, i2, XSIZE)
#define TALLOC_LL1(Z, SIZE) TALLOC_NEWINT(Sextra, Z, SIZE, j1, XXSIZE)
#define TALLOC_LL2(Z, SIZE) TALLOC_NEWINT(Sextra, Z, SIZE, j2, XXSIZE)
#define TALLOC_LLL1(Z, SIZE) TALLOC_NEWINT(Sextra, Z, SIZE,k1,XXXSIZE)
#define TALLOC_LLL2(Z, SIZE) TALLOC_NEWINT(Sextra, Z, SIZE,k2,XXXSIZE)

#define END_TALLOC_X1 END_TALLOC_NEW(a1)  //// // 116
#define END_TALLOC_X2 END_TALLOC_NEW(a2)////
#define END_TALLOC_X3 END_TALLOC_NEW(a3)////
#define END_TALLOC_XX1 END_TALLOC_NEW(b1) //// // 116
#define END_TALLOC_XX2 END_TALLOC_NEW(b2) ////
#define END_TALLOC_XXX1 END_TALLOC_NEW(c1) ////// 1116
#define END_TALLOC_XXX2 END_TALLOC_NEW(c2) //#


#define TALLOC_GATTER(z,SIZE) TALLOC_NEW(Sgatter, z, SIZE, z,XSIZE)
#define END_TALLOC_z END_TALLOC_NEW(z) // very crucial that z isn't freed!!
#define END_TALLOC_z1 END_TALLOC_NEW(z1)
#define END_TALLOC_z2 END_TALLOC_NEW(z2)


#define END_TALLOC_L1 END_TALLOC_NEW(i1)  //# # 16
//#define END_TALLOC_L2 END_TALLOC_NEW(i2)
#define END_TALLOC_LL1 END_TALLOC_NEW(j1)  // 116
//#define END_TALLOC_LL2 END_TALLOC_NEW(j2)
#define END_TALLOC_LLL1 END_TALLOC_NEW(k1)//# # 116
#define END_TALLOC_LLL2 END_TALLOC_NEW(k2)//#


int addShapeFct(model **Cov); /// ??

//typedef bool (*/* typusfct)(Types type); /\* h, cov, result *\/  */
/* //bool is_any(typusfct t, defn *C); */
/* //bool is_all(typusfct t, defn *C); */
/* bool isRObject(int type); */


// double gamma(double x);


// Formerly in <R_ext/Applic.h>LinkedTo: 

Rboolean fft_work(double *a, double *b, int nseg, int n, int nspn,
		  int isn, double *work, int *iwork,
		  int maxf, int kt, int m_fac, int *NFAC);/* TRUE: success // OK */
int fft_factor(int n, int *pmaxf, int *pmaxp, int *pkt, int *pm_fac, int*NFAC);




double searchInverse(covfct fct, model *cov, 
		     double start, double value, double releps);
double searchInverse(covfct fct, model *cov, 
		     double start, double min, double value, double releps);
// int searchFirstGreater(double *v, int len, double z);



#include "Machine.h"
#ifdef SCHLATHERS_MACHINE
#include "MachineSchlather.h"
#else 
#include "MachineOthers.h"
#endif

#ifdef RANDOMFIELDS_DEBUGGING
#include "MachineDebugging.h"
#endif


bool leading_spaces(model *cov, const char *character);
#define LPRINT if (!leading_spaces(cov, DOT)) {} else PRINTF
 


///////////////////////////////////////////////////////////////////////
// printing for debugging




void PSTOR(model *cov, gen_storage *x);
void PrintLoc(int level, location_type *loc, bool own);
void pmi(model *cov, char all, int level, int maxlevel); 
void pmi(model *cov, int maxlevel);
void pmiroot(model *cov, int maxlevel);
void pcl(int nr);
void pcl(model *cov);
void pcl();
void tree(model *cov, bool alle); // TREE
void printD(model *cov); //
void printI(model *cov); //




#define PMIL(Cov, ML) {				\
    PRINTF("\n(PMI '%s', line %d)", __FILE__, (int) __LINE__);	\
  pmi(Cov, ML);							\
}

#define PMI(Cov) PMIL(Cov, 999999)    
#define PMI0(Cov) PMIL(Cov, 0)
#define PMI1(Cov) PMIL(Cov, 1)
#define PMI2(Cov) PMIL(Cov, 2)
#define APMIL(Cov, ML) { PMIL(Cov, ML); assert(false); }	  
#define APMI(Cov)  APMIL(Cov, 999999)
#define APMI0(Cov) APMIL(Cov, 0)

#define PMIR(Cov) {						\
  PRINTF("\n(PMI '%s', line %d)", __FILE__, __LINE__);		\
  pmiroot(Cov, 999999);							\
}
#define PMIE(Cov) {					\
  PRINTF("\n\nPMIE '%s', line %d", __FILE__, __LINE__);	\
  PRINTF("\n%s level=%d err=%d (%s)\n\n", NAME(Cov), Cov->err_level, Cov->err, Cov->err_msg); \
  }
 
#define TREE(Cov) { PRINTF("\n(TREE '%s', line %d)\n", __FILE__, __LINE__); tree(Cov, true); }
#define TREE0(Cov) { PRINTF("\n(TREE '%s', line %d)\n", __FILE__, __LINE__); tree(Cov, false); }

    
#define PCL							\
    PRINTF("\nPCL '%s', line %d:  ", __FILE__, __LINE__);	\
    pcl
  

#define PLE PRINTF("\n(PLE '%s', line %d)", __FILE__, __LINE__); ple_
void ple_(model *cov);
void ple_(char *name);


void psys(system_type *sys);
void psys(model *cov);
#define PSYS(sys) { /* // */ 						\
	    PRINTF("\n("#sys" '%s', line %d)", __FILE__, __LINE__);	\
	    psys(sys);			       		\
	  }



extern char InternalName[],
  CovNames[MAXNRCOVFCTS][MAXCHAR],
  CovNickNames[MAXNRCOVFCTS][MAXCHAR];
extern int PL;
extern defn *DefList;
//extern int PrInL;				
extern int gaussmethod[Forbidden + 1];



// naechste 2 Zeilen nur notwendig, weil atan2 in Windows nicht
// ordentlich programmiert ist
//#ifdef WIN32
#define NEARBYWHAT 1e15
#define NEARBY(x)  (FLOOR((x) * (NEARBYWHAT) + 0.5) / (NEARBYWHAT))
//#else 
//#define NEARBY(x)  (x)
//#endif

#define EXTRA_Q if (cov->q == NULL) { QALLOC(4); cov->q[0] = cov->q[1] =cov->q[2] = cov->q[3] = RF_NAN;}

#define QVALUE cov->q[0]
#define QVALUE1 QVALUE
#define QVALUE2 cov->q[1]
#define QVALUE3 cov->q[2]
#define QVALUE4 cov->q[3]

#define SPLIT(I, MM, DIM, INDEX)			\
  int ii__ = I;						\
  for (int k__=0; k__<DIM##M1; k__++) {			\
    int j__ = ii__ / MM[k__];				\
    INDEX[k__] = ii__ % MM[k__];		      	\
    							\
    ii__ = j__;						\
  }							\
  INDEX[DIM##M1] = ii__

     // n=length(mm); k=rep(0,n); for (i in 0:(cumm[n]-1)) {ii = i; for (j in 1:(n-1)) {k[j] = ii %% mm[j]; ii = as.integer(ii / mm[j]) }; k[n] = ii; print(k)}


/// UNSORTIERT:

#define SET_NR(C,N) {							\
  set_nr(SYSOF(C), N);							\
  (C)->checked = false;							\
  (C)->IallowedDone = false;						\
  (C)->DallowedDone = false;						\
  (C)->initialised = false;						\
  (C)->zaehler = (C)->zaehler > 0 ? -(C)->zaehler : (C)->zaehler;	\
  /*also ? full_derivs rese_derivs   monotone logspeed finiterange loggiven hess
    taylorN tailN  ? */				\
}


bool allowedIfalse(model *cov);
bool allowedItrue(model *cov);
bool allowedDfalse(model *cov);
bool allowedDtrue(model *cov);
bool allowedD(model *cov);
bool allowedI(model *cov);
bool allowedIsubs(model *cov, model **sub, int z);

#define SET_CALLING(which, to){						\
    (which)->calling = to; /*if anything change, change also within CMbuild*/ \
    if ((to) != NULL) {							\
      (which)->root=(to)->root;						\
      (which)->base=(to)->base;						\
    }									\
}

#define SET_CALLING_NULL(which, c){		\
    (which)->calling = NULL;			\
    (which)->root=(c)->root;			\
    (which)->base=(c)->base;			\
  }

#define CHECKED {\
    assert(cov->checked);  /* da fehlt noch ein Befehlt?! */	\
}

#include "Error.h"


#define PLoffset -10

#if defined ZZ
#undef ZZ
#endif

#ifdef SCHLATHERS_MACHINE
#define ZZ printf("ZEILE %d in %s, %s\n",__LINE__, __FUNCTION__,__FILE__) //
#else
#define ZZ
#endif

bool parallel();

#define GETSTORAGE(Sx, From, Storage)			\
  assert((From) != NULL && (From)->S##Storage != NULL);	\
  Storage##_storage *Sx;				\
  Sx = (From)->S##Storage;				\
  assert(Sx != NULL)
#define getStorage(S, Storage)  GETSTORAGE(S, cov, Storage)

#define STOMODEL stoModel
#define GETSTOMODEL getStorage(STOMODEL, model)
#define DELSTOMODEL model_DELETE(&(cov->Smodel), cov);    
#define NEWSTOMODEL NEW_COV_STORAGE_SAVE(cov, model)
#define ONCE_NEWSTOMODEL ONCE_NEW_STORAGE(model)
#define COVMODELKEYS_GIVEN(cov) \
  ((cov)->Smodel != NULL && (cov)->Smodel->keys_given)
#define MODELKEYS_GIVEN COVMODELKEYS_GIVEN(cov) 
 

// if (L == NULL) ERR("register not initialised as likelihood method");

void Zero(model *cov, double *v);
void Zero(int *info, model *cov, double *v);

#define UNKNOWN_NUMBER_GRIDPTS -1

#endif

