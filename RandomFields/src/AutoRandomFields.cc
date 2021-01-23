#include "AutoRandomFields.h"

const char // muss in separater Zeile sein
// nachfolgende Zeile eingerueckt um 0 oder 2 Zeichen

*DOMAIN_NAMES[LAST_DOMAIN + 1] = { //RC
  "single variable", "kernel", "framework dependent", "submodel dependent",
  "framework or submodel dependent",
  "parameter dependent", "<keep copy>", "mismatch"},


  *OPTIMISER_NAMES[nOptimiser] = { // RC
  "optim", "optimx", "soma", "nloptr", // nloptr muss immer 4. sein
    "GenSA", "minqa", "pso", "DEoptim"},
  
  *NLOPTR_NAMES[nNLOPTR] = { // RC
    "NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L", 
    "NLOPT_GN_DIRECT_L_RAND", "NLOPT_GN_DIRECT_NOSCAL", 
    "NLOPT_GN_DIRECT_L_NOSCAL", "NLOPT_GN_DIRECT_L_RAND_NOSCAL", 
    "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L",
    "NLOPT_LN_PRAXIS", "NLOPT_GN_CRS2_LM",
    "NLOPT_LN_COBYLA", "NLOPT_LN_NELDERMEAD", 
    "NLOPT_LN_SBPLX", "NLOPT_LN_BOBYQA", 
    "NLOPT_GN_ISRES"},

  *LIKELIHOOD_NAMES[nLikelihood] = { // RC 
    "auto", "full", "composite", "tesselation"},

  *DUPLICATEDLOC_NAMES[nDuplicatedloc] = { // RC 
  "error", "risk error", "regard as repeated", "scatter"},

  *ISO_NAMES[LAST_ISO + 1] =  { // RC
    "isotropic", "space-isotropic", "vector-isotropic", "symmetric",
    "cartesian system", "gnomonic projection", // 5
    "orthographic projection",
    "spherical isotropic", "spherical symmetric", "spherical system", 
    "earth isotropic", // 10
    "earth symmetric",  "earth system",
    // "symmetric log-cartesian", "log-cartesian",
    "cylinder system",
    "non-dimension-reducing", "framework dependent", // 15
    "submodel dependent",
    "framework or submodel dependent",
    "parameter dependent",
    "<internal keep copy>", "<mismatch>"},
  
  *TYPE_NAMES[LASTTYPE + 1] = { // RC
    "tail correlation", "positive definite", "variogram", "negative definite",
    "point-shape function", "shape function", "trend", "distribution or shape",
    "of manifold type", "process", "method for Gauss process",  // 10
    "normed process (non-negative values with maximum value being 0 or 1)",
    "method for Brown-Resnick process",
    "Smith", "Schlather", "Poisson", "PoissonGauss", // 16
    "distribution family",
    "interface", 
    "mathematical operator", 
    "other type", // 20
    //
    "badly defined", "<same as previous>", "likelihood", "evaluation",// 24
    
    "mixed input", "string input", "<special I>", "<special II>",
    "<special III>", "<special IV>"//  31
  }, 

  *NEGATIVE_NAMES[-UNSET] = {
  "PARAM_DEP", "PREVMODEL_DEP", "SUBMODEL_DEP", "MISMATCH", "UNSET"
  },
  
  *MONOTONE_NAMES[MONOTONE_TOTAL] = { // RC
    "not set", "mismatch in monotonicity", "submodel dependent monotonicity",
    "previous model dependent monotonicity",
    "parameter dependent monotonicity",
    "not monotone", "monotone", "Gneiting-Schaback class", 
    "normal mixture", 
    "completely monotone",  
    "Bernstein"},

  *SORT_ORIGIN_NAMES[1 + (int) all_origins] = {
    "original model", "MLE conform", "all"
  },
 
  *MODE_NAMES[nr_modes] = {
    "careless", "sloppy", "easygoing", "normal", 
    "precise", "pedantic", "neurotic"},

  *OUTPUTMODE_NAMES[nr_output_modes] = {
    "sp", "RandomFields", "geoR"},

  *ANGLE_NAMES[last_angle_mode + 1] = {
    "radians", "degree"},

  *REPORTCOORD_NAMES[nr_reportcoord_modes] = {
    "always", "warn", "important", "never"},

  *UNITS_NAMES[nr_units] = {
    "", "km", "miles", "<time>", "<user defined>"},

  *COORD_SYS_NAMES[nr_coord_sys] = {
    "auto", "keep", "cartesian", "earth",
    "sphere", "gnomonic", "orthographic", "coordinate system mix"},
  
  *COORD_NAMES_GENERAL[2] = {
    "coords.x", "coords.T"},
  
  *CARTESIAN_SYS_NAMES[3] = {
    "cartesian", "gnomonic", "orthographic"},

  
  *TYPEOF_PARAM_NAMES[LASTSORTOF + 1] = {
    "variance", "covariance", "sd", "signed sd", "scale", 
    "diagonal", "aniso", "integer", "unspecified",  "trend",
    "nugget", "critical to estimate", "never verified",    
    "internally ignored", "forbidden to be estimated", "unkown",
    "variance (used internally)",
    "critical to estimate (used internally)", "only used by MLE", 
    "neither used by MLE nor returned",
    "standard", "include never returned", "internal", "all",
    "none"},

  *EQ_NAMES[nEQ_NAMES] = {"==", "!=", "<=", "<", ">=", ">"},
  
  *NAMES_OF_NAMES[nNamesOfNames] = {// see also LIST_OF_NAMES in init.general.cc
    "EQ_NAMES", "ISO_NAMES", // never change ordering !!, see also userinterface.cc  
    "DOMAIN_NAMES","TYPE_NAMES", 
    // from here on changings are ok
    "MONOTONE_NAMES",
    "MODE_NAMES", "OUTPUTMODE_NAMES", "REPORTCOORD_NAMES",
    "UNITS_NAMES", "COORD_SYS_NAMES", "CARTESIAN_SYS_NAMES",
    "TYPEOF_PARAM_NAMES"},
  
  *PROJECTION_NAMES[nProjections] = {
    "space", "time"},

  *RMCOV_X[nVAR2COV_METHODS] = { // c ==  COVARIATE_C_NAME 
  "origin", "center", "extremals", "all"},

  *FCTN_TYPE_NAMES[TOTAL_FCTN_TYPE] = {
  "Variogram", "Madogram", "Covariance", "Pseudomadogram", "Pseudovariogram"},

  *METHOD_NAMES[Forbidden+1]={"circulant", // RC
			     "cutoff",
			     "intrinsic",
			     "tbm", 
			     "spectral", //4
			     "direct",
			     "sequential",
			     "trend",
			     "average",
			     "nugget", //9
			     "coins",
			     "hyperplane",
			     "specific",
			     "any method", // nothing
			     "forbidden"},
  
  *POISSON_SCATTER_NAMES[nPOISSON_SCATTER] = {// must be in the ordering of
  "optimized", "standard", "any"}; // addPGS/local in extremes.cc and RMS.cc
 

  
