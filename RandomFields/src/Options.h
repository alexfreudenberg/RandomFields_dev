/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2017 -- 2017 Martin Schlather

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


#ifndef RFuser_H
#define RFuser_H 1


#define NAT_SCALE 1
#define MAX_CE_MEM 16777216

#define generalN 20
// IMPORTANT: all names of general must have at least 3 letters !!!
extern const char *general[generalN];
#define GENERAL_MODUS 0 
#define GENERAL_STORING 1
#define GENERAL_EXACTNESS 7
#define GENERAL_CLOSE 10
#define GENERAL_REPORTCOORD 16
struct general_param {
  char pch; /*  character shown after each simulation
    just for entertainment of the user
    except for "!", then the numbers are shown
	   */
  bool 
  na_rm_lines, vdim_close_together, storing,
  /* true: intermediate results are stored: might be rather memory consuming,
         but the simulation can (depending on the method chosen) be much faster
	 when done the second time with exactly the same parameters
	 do not forget to call DeleteAllKeys when all the simulation are done.
   false: intermediate results are not stored if SimulateRF is called or
         stored only until DoGauss is called.
	 minimum memory consumption, but might be rather slow if many
	 simulations are performed with exactly the same parameters
	 to be safe call DeleteAllKeys when all the simulations are done
	 Warning! init_foo may depend on GENERAL_STORING! This may cause problems
	 if GENERAL_STORING is changed between init_foo and do_foo from false to
	 true and if do_foo is called more than once! (SimulateRF is safe in
	 this case)

	 indifferent behaviour for the simulation methods if parameters are 
	 changed after each simulation.
	 In case, only MAXFIELDS different parameter sets are used, but after each 
	 step the parameter set is changed, use different keynr for each 
	 parametere set, and STORING==true, to get fast simulation 
  */
    sp_conform, /* should the simulation result be return in 
		   as an sp class or in the old form ?
		   --- getting obsolete in future --- todo
		*/
    detailed_output, returncall;
 
  int  mode, /* hightailing, fast, normal, save, pedantic */
    output, /* output mode, #alternative to mode that changes 
	       several other parameters;
	       vice versa, spConform sets 'output'
	     */
    reportcoord
  ; /* 
		   see convert.T PL_* and RF.h PL_*
		  */
  int naturalscaling;
  /* is called PracticalRange in R (see Chiles&Delfiner?!) 
   Note: RFparameters() allows only 0 and 1 as values!

   has an effect only if cov (and not only cov.local, e.g. Brownian Motion) is 
   defined
   0 : using the covariance function as defined
   1 : rescaling of cov fctn such that cov(1)~=0.05, if rescaling function 
       does not exist then failure, e.g. for Bessel model
   2 : exact or approximate value (e.g. expPLUScirc)
   3 : MLE (special needs taken into account, Long memory covariance functions
       have too Long tails, which are shortened (so threshold 0.05 is changed to
       higher values in this case))
  +10: if any of the above fails : numerical evaluation!
   else : using rescaling if rescaling function exists, otherwise without 
          rescaling
  */  

  int expected_number_simu,
    every,  // every 'every' line is announced if every>0

    // bool aniso;  // currerntly cannot be changed by the user !!

  // siehe InternalCov.cc line 350, simu.cc line 394
  /* expand scale immediately to aniso;
		 this allows higher flexibility in the choice of 
		 simulation methods (tbm2 does not like dimension reduction),
		 but it is slower
	      */
    Ttriple, 
    set, seed_incr, seed_sub_incr, duplicated_loc;
 
  double gridtolerance;
  usr_bool exactness;

};

#define general_START \
  {pch[NM],								\
      false, false, false /*storing*/, true, /* 5 */ \
      false, false, 						\
      startmode/* mode */ , output_sp, reportcoord_warnings,		\
      false /* naturalscaling */,					\
      1, 0,								\
      NA_INTEGER, 0, 0,	101101,	DUPLICATEDLOC_RISKERROR,		\
      1e-6, exactness[NM]						\
     }


#define gaussN 6
extern const char *gauss[gaussN];
#define GAUSS_BEST_DIRECT 3
#define GAUSS_BOXCOX_OPTION 5
struct gauss_param{
  usr_bool stationary_only; // logical + NA
  double
    approx_zero, 
    boxcox[2 * MAXBOXCOXVDIM];
  bool paired, loggauss;
  int direct_bestvariables;
};
#define gauss_START  {Nan, 0.05,					\
      {RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0,		\
	  RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0, RF_INF, 0},	\
      false, false, 1200}

#define MINSPLITN(locmaxn) ((locmaxn) / 40)
#define MEDSPLITN(locmaxn) ((locmaxn) / 8)
#define MAXSPLITN(locmaxn) (((locmaxn) * 5) / 8)

#define krigeN 5
#define KRIGE_SPLITN 2
extern const char *krige[krigeN];
struct krige_param {
  bool ret_variance,  fillall
     ;
  int locmaxn,
    locsplitn[3], // 0:when splitting is done; 1:min pts in a neighbourhood ; 2:max pts when still neighbourhoods are included
    locsplitfactor;
};
#define krige_START {false, true,					\
      locmaxn[NM],							\
      {MINSPLITN(locmaxn[NM]), MEDSPLITN(locmaxn[NM]), MAXSPLITN(locmaxn[NM])},\
      2}


#define TRIVIALSTRATEGY 0
#define STRATEGY_PARTIAL 1
#define LASTSTRATEGY STRATEGY_PARTIAL


#define CEN 12
extern const char *CE[CEN];
struct ce_param {
  bool force, useprimes, dependent;
  char strategy;
  int trials, maxgridsize, maxmem;
  double maxGB,
     tol_re, tol_im, mmin[MAXCEDIM],
    approx_grid_step;
};
#define ce_START				\
  {ce_force[NM], true, ce_dependent[NM],	\
      TRIVIALSTRATEGY,				\
      3, MAX_CE_MEM, MAXINT,						\
      1, ce_tolRe[NM], ce_tolIm[NM], {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, \
      ce_approx_step[NM]}

#define spectralN 4
#define SPECTRAL_PROPFACTOR 2
extern const char * spectral[spectralN];
struct spectral_param {
  bool grid;
  double prop_factor, sigma;
  int lines[MAXTBMSPDIM];
};
#define spectral_START {sp_grid[NM], 50, 0.0, {2500, 2500, 2500, 2500}}



#define pTBMN 9
extern const char * pTBM[pTBMN];
struct tbm_param {
  bool grid;
  int tbmdim, fulldim, points,
    lines[MAXTBMSPDIM];          // number of lines simulated
  usr_bool layers; // 0, 1, NA
  double linesimufactor, /*factor by which simulation on the line is finer 
			     than the grid */
    linesimustep,      /* grid lag defined absolutely */
    center[MAXTBMSPDIM];
};
#define tbm_START \
  {sp_grid[NM], -2, 3, 0, {1, 60, 500}, Nan, tbm_linefactor[NM], 0.0, \
   {RF_NA, RF_NA, RF_NA, RF_NA}}


#define MAX_DIRECT_MAXVAR 30000
#define DIRECT_ORIG_MAXVAR 12000

#define directN 1
#define DIRECT_MAXVAR_PARAM 0
extern const char *direct[directN];
struct direct_param {
  int maxvariables;
};
#define direct_START  { DIRECT_ORIG_MAXVAR }


#define sequN 2
extern const char * sequ[sequN];
struct sequ_param{
  int back, initial;
};
#define sequ_START {10, -10}



struct ave_param {
};
#define ave_START {}


#define pnuggetN 1
extern const char * pnugget[pnuggetN];
struct nugget_param {
  double tol;
};
#define nugget_START {nugget_tol[NM]}


#define mppN 6
extern const char * mpp[mppN];
struct mpp_param {
  int  n_estim_E,
    scatter_max[MAXMPPDIM];
  double intensity[MAXMPPDIM], // intensity factor for e.g. unif_initu
    about_zero,
    shape_power,
    scatter_step[MAXMPPDIM];
};
#define mpp_START					\
  {50000,						\
   {NA_INTEGER, NA_INTEGER, NA_INTEGER, NA_INTEGER},			\
   {mpp_intensity[NM], mpp_intensity[NM], mpp_intensity[NM],mpp_intensity[NM]},\
     mpp_zero[NM], /* about zero */					\
       2.0,								\
	 {RF_NA, RF_NA, RF_NA, RF_NA}					\
  }



#define HYPER_UNIFORM  0
#define HYPER_FRECHET  1 
#define HYPER_BERNOULLI 2

#define hyperN 4
extern const char * hyper[hyperN];
struct hyper_param {
  int superpos, maxlines, mar_distr;
  double mar_param;
};
#define hyper_START {700, 1000, HYPER_UNIFORM, RF_NA}


#define extremeN 12
extern const char * extreme[extremeN];
#define EXTREME_FLAT 5
struct extremes_param {
  usr_bool flathull;
  int maxpoints, check_every, min_n_zhou, max_n_zhou, mcmc_zhou, scatter_method;
  double standardmax,
    GEV_xi, density_ratio, eps_zhou, min_shape_gumbel;
};
#define extreme_START				\
  {False, MAXINT, 30,			\
      1000, 10000000, 20, POISSON_SCATTER_ANY ,	\
      max_max_gauss[NM],			\
      1.0, 0.0, 0.01, -1e15 /* == -inf */	\
}


#define brN 7
extern const char * br[brN];
struct br_param {
  int BRmaxmem, BRvertnumber, BRoptim, deltaAM;
  double BRmeshsize, BRoptimtol, variobound;
};
#define br_START	       \
  {10000000, 7, 1, 300, \
      0.1, 0.01, 8.0}     \


#define distrN 9
extern const char * distr[distrN];
struct distr_param{
  double safety, minsteplen, innermin, outermax;
  int maxsteps, parts, maxit, mcmc_n, repetitions;
};
#define distr_START {0.08, 0, 1e-20, 1e5,				\
      DISTMAXSTEPS, 8, 20, 15, 1000}  // distr (rectangular) // todo should be 500 and better algorithm for approximation!



#define MINCLIQUE(locmaxn) ((locmaxn) / 25)
#define MEDCLIQUE(locmaxn) ((locmaxn) / 5)
#define MAXCLIQUE(locmaxn) (((locmaxn) * 3) / 5)

#define FLAT_UNDETERMINED -1
#define fitN 43
#define FIT_BC_LB 10
#define FIT_BC_UB 11
#define FIT_MAXNEIGHBOUR 19
#define FIT_CLIQUE 20
extern const char * fit[fitN];
struct fit_param{
  double bin_dist_factor, upperbound_scale_factor, lowerbound_scale_factor, 
    lowerbound_scale_LS_factor, upperbound_var_factor, lowerbound_var_factor, 
  //lowerbound_sill, 
    scale_max_relative_factor, minbounddistance, 
    minboundreldist,  
    BC_lambdaLB[2 * MAXBOXCOXVDIM],
    BC_lambdaUB[2 * MAXBOXCOXVDIM],
    scale_ratio,
    min_diag, pgtol, pgtol_recall,
    factr, factr_recall, emp_alpha ;// logical + NA // 19
  int approximate_functioncalls,
     critical, 
    n_crit, locmaxn, 
    locsplitn[3], // 0:when splitting is done; 1:min pts in a neighbourhood ; 2:max pts when still neighbourhoods are included
    locsplitfactor, smalldataset, optimiser, suboptimiser,
    algorithm, likelihood, split, addNAlintrend, trace; // 
  usr_bool estimate_variance;
  bool  use_naturalscaling, onlyuser, 
   reoptimise, ratiotest_approx,
    split_refined, cross_refit, suggesting_bounds, return_hessian,
    standardizedL;      // 8
  char lengthshortname; // 1..255
};
#define fit_START\
  0.4, 3.0, 5.0, 3.0, 10.0, 1e4, /* 6 */				\
    1000.0, 0.001, 0.02,  /* 11 */					\
  {RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF,	\
      RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, \
      RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, RF_NEGINF, \
      RF_NEGINF, RF_NEGINF},						\
	{RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, \
	    RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, RF_INF, \
	    RF_INF, RF_INF, RF_INF, RF_INF},				\
    0.1, 1e-7, fit_pgtol[NM], fit_pgtol_recall[NM], /* */		\
    fit_factr[NM], fit_factr_recall[NM], PSEUDOVARIOGRAM,		\
    /* int: */								\
    50,  0 /* critical */,            /* 6 */				\
    5 /* ncrit */ , maxclique[NM],					\
	    {MINCLIQUE(maxclique[NM]), MEDCLIQUE(maxclique[NM]),	\
		MAXCLIQUE(maxclique[NM]) }, 2, 2000 /* small */,  /* 11 */ \
    0, 0, 0 /* UNSET algorithm */, 0, fit_split[NM], 1, 0, /*  */	\
    /* usr_bool */ Nan,							\
	  /* bool */ false, false, fit_reoptimise[NM],  /* 5 */ \
    fit_ratiotest_approx[NM], true, fit_cross_refit[NM], false, true, false, \
    12 // fit

#define empvarioN 9
extern const char * empvario[empvarioN];
struct empvario_param{
  double phi0, theta0, tol,
    bins[MAXBINS], deltaT[2];
  
  bool fft, halveangles;
  int nbins, nphi, ntheta;
};
#define empvario_START							\
  {0.0, 0.0, 1e-13,							\
      { 20 }, {0, 0},							\
	     false, true,						\
	       1,  0, 0}

#define guiN 3
extern const char * gui[guiN];
#define GUI_SIZE 2
struct gui_param{
  bool alwaysSimulate;
  int method, size[2];
};
#define gui_START {true, CircEmbed, {1024, 64}}

#define graphicsN 16
extern const char *graphics[graphicsN];
#define GRAPHICS_UPTO 3
struct graphics_param{
  usr_bool always_open, always_close;
  bool split_screen, close_screen, onefile, grDefault;
  double height, width, resolution, xlim;
  int PL, increase_upto[2], number, maxchar, npoints[2];
  char filename[100];
};

#define graphics_START  { True, Nan,		\
      true, true, false, false,			\
      6.0, RF_NA, 72.0,	3.0,			\
      0, {3, 4}, 0, 20, {100, 200}, ""}


#define registersN 3
extern const char *registers[registersN];
struct registers_param {
  int keynr, predict, likelihood;
};
#define register_START {0, MODEL_PREDICT, MODEL_USER}


#define internalN 2
extern const char * internals[internalN];
#define INTERNALS_DO_TESTS 0
#define INTERNALS_EX_RED 1
struct internal_param{ 
  // if changed, CHANGE ALSO RestWarnings in 'userinterfaces.cc';
  bool do_tests;
  int examples_reduced;
};
#define internal_START			\
  {DO_TESTS, 0}


#define NO_NOTE 0
#define NOTE_ONCE 1
#define NOTE_WITHOUT_HINT 2
#define NOTE_WITH_HINT 3

#define messagesN 28
extern const char * messages[messagesN];
#define MESSAGES_NEWANISO 2
#define MESSAGES_ONGRID 8
#define MESSAGES_COORD_CHANGE 11
#define MESSAGES_ZENIT 13
#define MESSAGES_RAW 23
struct messages_param{ 
  // if changed, CHANGE ALSO RestWarnings in 'userinterfaces.cc';
  bool
  warn_oldstyle, help_newstyle, warn_Aniso, note_ambiguous, help_normal_mode,
    warn_mode, warn_scale, warn_on_grid, warn_ambiguous, 
    note_no_fit, note_aspect_ratio, warn_coord_change, help_color_palette,
    warn_zenit, // not used anymore
    warn_constant, warn_negvar,  help_onlyvar, 
    warn_modus_operandi, warn_singlevariab, help_mle,
    warn_raw_covariates, help_addNA, help_help;
  usr_bool warn_mathdef  // obsolete?!
    ;//   
  int  note_detection, note_coordinates, warn_seed, vec_len;
};
#define messages_START			\
  {true, true, true, false, true,		\
      true, true, true, true,		\
      true, true, true, true, true,		\
      true, true, true, 		\
      true, true, true, \
      true, true, true,				\
      Nan,					\
      NOTE_WITH_HINT, NOTE_WITH_HINT, NOTE_WITH_HINT, 10}


#define coordsN 20
#define COORDS_XYZNOTATION 0
#define COORDS_DATAIDX 5
#define COORDS_DATANAMES 6
#define COORDS_XIDX 7
#define COORDS_XNAMES 8
#define ZENIT 10
extern const char *coords[coordsN];
struct coords_param {
  usr_bool xyz_notation; // 0; + Exception for RFcov
  double zenit[2]; // 8
  coord_sys_enum coord_system, // 1
    new_coord_system; // 7
  char
     newunits[MAXCOORDNAMES][MAXUNITSCHAR], // 2; only to read
    curunits[MAXCOORDNAMES][MAXUNITSCHAR], // 3
    varunits[MAXDATANAMES][MAXUNITSCHAR], // 4

  // 2 user variables for data.frames
    data_names[MAXDATANAMES][MAXCHAR],// (5C)
    x_names[MAXCOORDNAMES][MAXCHAR], // (6 C)
    data_col_initial[MAXCHAR],
    coord_initial[MAXCHAR],
    cartesian_names[4][MAXCHAR],
    earthcoord_names[4][MAXCHAR];
  // auxiliary variables 
  int
    data_nr_names, //  5; needed by data_names and x_names
    x_nr_names, //  5; needed by data_names and x_names
    data_idx[2],
    x_idx[2], // (5B, 6B) alternatives for data_names and x_namwe
    max_col, max_coord;
  bool polar_coord, // 9
    allow_earth2cart; // 11
  angle_modes anglemode; // 10
};
#define coords_START 							\
  Nan, {1, RF_NA}, coord_auto, coord_keep,				\
    {""}, {""}, {""}, {""}, {""},				\
    "X", "V", {"x", "y", "z", "T"}, {"longitude","latitude","height","time"},	\
    0, 0, {NA_INTEGER, NA_INTEGER}, {NA_INTEGER, NA_INTEGER}, 99, 4,	\
    false, false, radians


#define specialN 1
extern const char * special[specialN];
struct special_param {
   int multcopies;
};
#define special_START {20}


#define obsoleteN 10
extern const char * obsolete[obsoleteN];


#define ignoreN 1
// trend

#define prefixN 25
struct globalparam{
  general_param general;
  gauss_param gauss;
  krige_param krige;
  ce_param ce;
  spectral_param spectral;


  tbm_param tbm;
  direct_param direct;
  sequ_param sequ;
  ave_param ave;
  nugget_param nugget;
  mpp_param mpp;
   hyper_param hyper;
  extremes_param extreme;
  br_param br;
  distr_param distr;
  fit_param fit;
 
  empvario_param empvario;
  gui_param gui;
  graphics_param graphics;
  registers_param registers;
  internal_param internal;
  messages_param messages;
  coords_param coords;
  special_param special;
};

extern const char * prefixlist[prefixN], **all[prefixN];
extern int allN[prefixN];

extern globalparam GLOBAL;

#endif
