
#ifndef rfutils_options_H
#define rfutils_options_H 1

#include <R.h>
#include <Rdefines.h>
#include "RFU.h"
#include "AutoRandomFieldsUtilsLocal.h"

#if defined __WIN32__ || defined (__APPLE__) || defined(__sun)
#else
#define TIME_AVAILABLE 1
#endif

#define basicN 12
// IMPORTANT: all names of basic must be have least 3 letters !!!
extern const char *basic[basicN];
#define BASIC_WARN_OPTION 9
typedef // benoetigt
struct basic_param {
  int  
  Rprintlevel, Cprintlevel, seed, cores, warn_unknown_option,
    LaMaxTakeOwn;
  bool skipchecks, asList /* hidden:verbose */, kahanCorrection, helpinfo,
   warn_parallel;
  la_modes la_usr, la_mode;
} basic_param;
#define basic_START \
  { R_PRINTLEVEL, C_PRINTLEVEL, NA_INTEGER, INITCORES, WARN_UNKNOWN_OPTION_ALL,\
      MAXINT, 							\
      false, true, false, true, true,					\
      LA_AUTO, LA_R							\
      }

extern const char * InversionNames[nr_InversionMethods];

#define SOLVE_SVD_TOL 3
#define solveN 20
typedef // benoetigt
struct solve_param {
  usr_bool sparse, pivot_check;
  bool det_as_log, pivot_partialdet;
  double spam_tol, spam_min_p[2], svd_tol, eigen2zero, pivot_relerror,
    max_deviation, max_reldeviation;
  InversionMethod Methods[SOLVE_METHODS];
  int spam_min_n[2], spam_sample_n, spam_factor, pivotsparse, max_chol,
    max_svd,
     actual_size,
    *pivot_idx, pivot_idx_n,//permutation; phys+logi laenge
    tinysize;
  //  bool tmp_delete;
  pivot_modes actual_pivot,pivot_mode;
} solve_param;
#ifdef SCHLATHERS_MACHINE
#define svd_tol_start 1e-08
#else
#define svd_tol_start 0
#endif
#define solve_START							\
  False, False, true, false,						\
    2.220446e-16, {0.8, 0.9}, svd_tol_start, 1e-12, 1e-11,		\
    1e-10, 1e-10,							\
    {NoInversionMethod,  NoFurtherInversionMethod},			\
    {400, 10000}, 500, 4294967, PIVOTSPARSE_MMD, 16384,			\
    10000,  /* never change -- see RFoptions.Rd */			\
    0, NULL, 0, 3,							\
    PIVOT_UNDEFINED, PIVOT_NONE
    
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
//int own_chol_up_to();
void SetLaMode();
void SetLaMode(la_modes);


#endif
