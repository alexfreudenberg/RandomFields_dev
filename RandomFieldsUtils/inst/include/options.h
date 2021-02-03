
#ifndef rfutils_options_H
#define rfutils_options_H 1

#include <R.h>
#include <Rdefines.h>
#include "RFU.h"


#define basicN 12
// IMPORTANT: all names of basic must be have least 3 letters !!!
extern const char *basic[basicN];
#define BASIC_WARN_OPTION 9
#define BASIC_USEGPU 10
typedef // benoetigt
struct basic_param {
  int  
  Rprintlevel, Cprintlevel, seed, cores, warn_unknown_option;
  bool skipchecks, asList /* hidden:verbose */, kahanCorrection, helpinfo,
    useGPU, warn_parallel;
} basic_param;
#define basic_START \
  { R_PRINTLEVEL, C_PRINTLEVEL, NA_INTEGER, INITCORES, WARN_UNKNOWN_OPTION_ALL,\
      false, true, false, true,  GPUavailable, true			\
      }


extern const char * InversionNames[nr_InversionMethods];

#define SOLVE_SVD_TOL 3
#define solveN 20
typedef // benoetigt
struct solve_param {
  usr_bool sparse, pivot_check;
  bool det_as_log;
  double spam_tol, spam_min_p[2], svd_tol, eigen2zero, pivot_relerror,
    max_deviation, max_reldeviation;
  InversionMethod Methods[SOLVE_METHODS];
  int spam_min_n[2], spam_sample_n, spam_factor, pivotsparse, max_chol,
    max_svd, pivot, actual_pivot, actual_size,
    *pivot_idx, pivot_idx_n,//permutation; phys+logi laenge
    pivotMaxTakeOwn;
  //  bool tmp_delete;
} solve_param;
#ifdef SCHLATHERS_MACHINE
#define svd_tol_start 1e-08
#else
#define svd_tol_start 0
#endif
#define solve_START							\
  { False, False, true,							\
      DBL_EPSILON, {0.8, 0.9}, svd_tol_start, 1e-12, 1e-11,		\
    1e-10, 1e-10,							\
  {NoInversionMethod,  NoFurtherInversionMethod},			\
    {400, 10000}, 500, 4294967, PIVOTSPARSE_MMD, 16384,			\
	10000, PIVOT_NONE, /* never change -- see RFoptions.Rd */	\
		    PIVOT_UNDEFINED, 0, NULL, 0, NA_INTEGER}
extern const char * solve[solveN];

typedef // benoetigt
struct utilsparam{
  basic_param basic;
  solve_param solve;
} utilsparam;


#define LEN_OPTIONNAME 201 // zwingend ungerade
typedef void (*setoptions_fctn) (int, int, SEXP, char[LEN_OPTIONNAME],
				 bool, bool);
typedef void (*getoptions_fctn) (SEXP, int, bool);
typedef void (*finalsetoptions_fctn) ();
typedef void (*deleteoptions_fctn) (bool);


#define R_PRINTLEVEL 1
#define C_PRINTLEVEL 1
#ifdef SCHLATHERS_MACHINE
#define INITCORES 4
#else
#define INITCORES 1
#endif


#define ADD(ELT) SET_VECTOR_ELT(sublist, k++, ELT)
#define ADDCHAR(ELT) x[0] = ELT; ADD(ScalarString(mkChar(x)))


//int own_chol_up_to(int size, int maxtime);
int own_chol_up_to();


#endif
