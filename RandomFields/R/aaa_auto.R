# This file has been created automatically by 'rfGenerateConstants'


 ## from  ../../RandomFieldsUtils/RandomFieldsUtils/src/AutoRandomFieldsUtils.h

 MAXUNITS 	<- as.integer(4)
 MAXCHAR 	<- as.integer(18)
 RFOPTIONS 	<- "RFoptions"

 WARN_UNKNOWN_OPTION_ALL 	<- as.integer(3)
 WARN_UNKNOWN_OPTION_SINGLE 	<- as.integer(2)
 WARN_UNKNOWN_OPTION_CAPITAL 	<- as.integer(1)
 WARN_UNKNOWN_OPTION_NONE 	<- as.integer(0)
 WARN_UNKNOWN_OPTION 	<- as.integer(10000)
 WARN_UNKNOWN_OPTION_CONDSINGLE 	<- as.integer((WARN_UNKNOWN_OPTION_SINGLE-WARN_UNKNOWN_OPTION))
 WARN_UNKNOWN_OPTION_DEFAULT 	<- as.integer(WARN_UNKNOWN_OPTION_ALL)



 ## from  src/AutoRandomFields.h

 AUTHOR 	<- "martin.schlather@math.uni-mannheim.de.\\"
 CONTACT 	<- " Please contact the maintainer martin.schlather@math.uni-mannheim.de.\n"

 MAXCHAR_RF 	<- as.integer(MAXCHAR)
 METHODMAXCHAR 	<- as.integer(MAXCHAR)
 MAXPARAM 	<- as.integer(20)
 MAXSUB 	<- as.integer(10)

 MAXCEDIM 	<- as.integer(13)
 MAXTBMSPDIM 	<- as.integer(4)
 MAXMPPDIM 	<- as.integer(4)

 MAXBOXCOXVDIM 	<- as.integer(10)

 MAXHYPERDIM 	<- as.integer(4)
 MAXVARIODIM 	<- as.integer(20)
 MAXACCEPTED 	<- as.double(1e300)
 MAXVARIANTS 	<- as.integer(6)
 MAXSYSTEMS 	<- as.integer(1)
 MAXDATANAMES 	<- as.integer(5)

 MAXMPPVDIM 	<- as.integer(9)

 MAXCOORDNAMES 	<- as.integer(4)
 units_none 	<- as.integer(0)
 units_km 	<- as.integer(1)
 units_miles 	<- as.integer(2)
 units_time 	<- as.integer(3)
 units_user 	<- as.integer(4)

 nr_units 	<- as.integer((units_user+1))

 coord_auto 	<- as.integer(0)
 coord_keep 	<- as.integer(1)
 cartesian 	<- as.integer(2)
 earth 	<- as.integer(3)
 sphere 	<- as.integer(4)
 gnomonic 	<- as.integer(5)
 orthographic 	<- as.integer(6)
 coord_mix 	<- as.integer(7)


 nr_coord_sys 	<- as.integer((coord_mix+1))

 reportcoord_always 	<- as.integer(0)
 reportcoord_warnings_orally 	<- as.integer(1)
 reportcoord_warnings 	<- as.integer(2)
 reportcoord_none 	<- as.integer(3)

 nr_reportcoord_modes 	<- as.integer((reportcoord_none+1))

 careless 	<- as.integer(0)
 sloppy 	<- as.integer(1)
 easygoing 	<- as.integer(2)
 normal 	<- as.integer(3)
 precise 	<- as.integer(4)
 pedantic 	<- as.integer(5)
 neurotic 	<- as.integer(6)

 nr_modes 	<- as.integer((neurotic+1))

 output_sp 	<- as.integer(0)
 output_rf 	<- as.integer(1)
 output_geor 	<- as.integer(2)

 nr_output_modes 	<- as.integer((output_geor+1))

 radians 	<- as.integer(0)
 degree 	<- as.integer(1)

 last_angle_mode 	<- as.integer(degree)

 noConcerns 	<- as.integer(0)
 neverRaw 	<- as.integer(1)
 neverGlobalXT 	<- as.integer(2)
 ignoreValues 	<- as.integer(3)
 unsetConcerns 	<- as.integer(4)


 PARAM_DEP 	<- as.integer(-1)
 PREVMODEL_DEP 	<- as.integer(-2)
 SUBMODEL_DEP 	<- as.integer(-3)
 MISMATCH 	<- as.integer(-4)
 UNSET 	<- as.integer(-5)

 XONLY 	<- as.integer(0)
 KERNEL 	<- as.integer(1)
 PREVMODEL_D 	<- as.integer(2)
 SUBMODEL_D 	<- as.integer(3)
 PREV_SUB_D 	<- as.integer(4)
 PARAMDEP_D 	<- as.integer(5)
 KEEPCOPY_DOM 	<- as.integer(6)
 DOMAIN_MISMATCH 	<- as.integer(7)


 FIRST_DOMAIN 	<- as.integer(XONLY)
 LAST_DOMAINUSER 	<- as.integer(KERNEL)
 LAST_DOMAIN 	<- as.integer(DOMAIN_MISMATCH)

RC_ISOTROPIC <- ISOTROPIC 	<- as.integer(0)
RC_DOUBLEISOTROPIC <- DOUBLEISOTROPIC 	<- as.integer(1)
 VECTORISOTROPIC 	<- as.integer(2)
 SYMMETRIC 	<- as.integer(3)
RC_CARTESIAN_COORD <- CARTESIAN_COORD 	<- as.integer(4)
RC_GNOMONIC_PROJ <- GNOMONIC_PROJ 	<- as.integer(5)
RC_ORTHOGRAPHIC_PROJ <- ORTHOGRAPHIC_PROJ 	<- as.integer(6)
 SPHERICAL_ISOTROPIC 	<- as.integer(7)
 SPHERICAL_SYMMETRIC 	<- as.integer(8)
 SPHERICAL_COORD 	<- as.integer(9)
RC_EARTH_ISOTROPIC <- EARTH_ISOTROPIC 	<- as.integer(10)
 EARTH_SYMMETRIC 	<- as.integer(11)
 EARTH_COORD 	<- as.integer(12)
 CYLINDER_COORD 	<- as.integer(13)
RC_UNREDUCED <- UNREDUCED 	<- as.integer(14)
 PREVMODEL_I 	<- as.integer(15)
 SUBMODEL_I 	<- as.integer(16)
 PREV_SUB_I 	<- as.integer(17)
 PARAMDEP_I 	<- as.integer(18)
 KEEPCOPY_ISO 	<- as.integer(19)
 ISO_MISMATCH 	<- as.integer(20)


 FIRST_CARTESIAN 	<- as.integer(ISOTROPIC)
 FIRST_ISOUSER 	<- as.integer(FIRST_CARTESIAN)
 LAST_REDUCEDXDIM_CART 	<- as.integer(DOUBLEISOTROPIC)
 FIRST_PROJECTION 	<- as.integer(GNOMONIC_PROJ)
 LAST_PROJECTION 	<- as.integer(ORTHOGRAPHIC_PROJ)
 LAST_CARTESIAN 	<- as.integer(LAST_PROJECTION)
 FIRST_SPHERICAL 	<- as.integer(SPHERICAL_ISOTROPIC)
 FIRST_ANYSPHERICAL 	<- as.integer(FIRST_SPHERICAL)
 LAST_TRUESPHERICAL 	<- as.integer(SPHERICAL_COORD)
 FIRST_EARTH 	<- as.integer(EARTH_ISOTROPIC)
 LAST_EARTH 	<- as.integer(EARTH_COORD)
 LAST_ANYSPHERICAL 	<- as.integer(LAST_EARTH)

 LAST_ISOUSER 	<- as.integer(UNREDUCED)
 LAST_ISO 	<- as.integer(ISO_MISMATCH)

 MON_UNSET 	<- as.integer(UNSET)
 MON_MISMATCH 	<- as.integer(MISMATCH)
 MON_SUB_DEP 	<- as.integer(SUBMODEL_DEP)
 MON_PREV_DEP 	<- as.integer(PREVMODEL_DEP)
 MON_PARAMETER 	<- as.integer(PARAM_DEP)
 NOT_MONOTONE 	<- as.integer(0)
 MONOTONE 	<- as.integer(1)
 GNEITING_MON 	<- as.integer(2)
 NORMAL_MIXTURE 	<- as.integer(3)
 COMPLETELY_MON 	<- as.integer(4)
 BERNSTEIN 	<- as.integer(5)


 MONOTONE_TOTAL 	<- as.integer((BERNSTEIN-MON_UNSET+1))

 MAXFIELDS 	<- as.integer(10)
 MODEL_USER 	<- as.integer((MAXFIELDS+0))
 MODEL_COV 	<- as.integer((MAXFIELDS+1))
 MODEL_COVMATRIX 	<- as.integer((MAXFIELDS+2))
 MODEL_VARIOGRAM 	<- as.integer((MAXFIELDS+3))
 MODEL_PSEUDO 	<- as.integer((MAXFIELDS+4))
 MODEL_FCTN 	<- as.integer((MAXFIELDS+5))
 MODEL_DISTR 	<- as.integer((MAXFIELDS+6))
 MODEL_CALC 	<- as.integer((MAXFIELDS+7))

 LAST_MODEL_USER 	<- as.integer((MAXFIELDS+9))
 FIRST_INTERNAL 	<- as.integer((LAST_MODEL_USER+1))
 MODEL_AUX 	<- as.integer((FIRST_INTERNAL+0))
 MODEL_INTERN 	<- as.integer((FIRST_INTERNAL+1))
 MODEL_SPLIT 	<- as.integer((FIRST_INTERNAL+2))
 MODEL_GUI 	<- as.integer((FIRST_INTERNAL+3))
 MODEL_MLE 	<- as.integer((FIRST_INTERNAL+4))
 MODEL_MLESPLIT 	<- as.integer((FIRST_INTERNAL+5))
 MODEL_LSQ 	<- as.integer((FIRST_INTERNAL+6))
 MODEL_BOUNDS 	<- as.integer((FIRST_INTERNAL+7))
 MODEL_PREDICT 	<- as.integer((FIRST_INTERNAL+8))
 MODEL_ERR 	<- as.integer((FIRST_INTERNAL+10))
 MODEL_MAX 	<- as.integer(MODEL_ERR)

 TcfType 	<- as.integer(0)
 PosDefType 	<- as.integer(1)
 VariogramType 	<- as.integer(2)
 NegDefType 	<- as.integer(3)
 PointShapeType 	<- as.integer(4)
 ShapeType 	<- as.integer(5)
 TrendType 	<- as.integer(6)
 RandomOrShapeType 	<- as.integer(7)
 ManifoldType 	<- as.integer(8)
 ProcessType 	<- as.integer(9)
 GaussMethodType 	<- as.integer(10)
 NormedProcessType 	<- as.integer(11)
 BrMethodType 	<- as.integer(12)
 SmithType 	<- as.integer(13)
 SchlatherType 	<- as.integer(14)
 PoissonType 	<- as.integer(15)
 PoissonGaussType 	<- as.integer(16)
 RandomType 	<- as.integer(17)
 InterfaceType 	<- as.integer(18)
 MathDefType 	<- as.integer(19)
 OtherType 	<- as.integer(20)
 BadType 	<- as.integer(21)
 SameAsPrevType 	<- as.integer(22)
 LikelihoodType 	<- as.integer(23)
 EvaluationType 	<- as.integer(24)
 MixedInputType 	<- as.integer(25)
 CharInputType 	<- as.integer(26)
 NN1 	<- as.integer(27)
 NN2 	<- as.integer(28)
 NN3 	<- as.integer(29)
 NN4 	<- as.integer(30)


 LASTNN 	<- as.integer(NN4)
 LASTTYPE 	<- as.integer(NN4)

 nOptimiser 	<- as.integer(9)
 nNLOPTR 	<- as.integer(15)
 nLikelihood 	<- as.integer(4)
 nDuplicatedloc 	<- as.integer(4)
 nNamesOfNames 	<- as.integer(12)
 nVAR2COV_METHODS 	<- as.integer(4)
 nProjections 	<- as.integer(2)

 DUPLICATEDLOC_ERROR 	<- as.integer(0)
 DUPLICATEDLOC_RISKERROR 	<- as.integer(1)
 DUPLICATEDLOC_REPEATED 	<- as.integer(2)
 DUPLICATEDLOC_SCATTER 	<- as.integer(3)
 DUPLICATEDLOC_LASTERROR 	<- as.integer(DUPLICATEDLOC_RISKERROR)

 DetTrendEffect 	<- as.integer(0)
 FixedTrendEffect 	<- as.integer(1)
 FixedEffect 	<- as.integer(2)
 DataEffect 	<- as.integer(3)
 RandomEffect 	<- as.integer(4)
 ErrorEffect 	<- as.integer(5)
 effect_error 	<- as.integer(6)


 VARPARAM 	<- as.integer(0)
 SIGNEDVARPARAM 	<- as.integer(1)
 SDPARAM 	<- as.integer(2)
 SIGNEDSDPARAM 	<- as.integer(3)
 SCALEPARAM 	<- as.integer(4)
 DIAGPARAM 	<- as.integer(5)
 ANISOPARAM 	<- as.integer(6)
 INTEGERPARAM 	<- as.integer(7)
 ANYPARAM 	<- as.integer(8)
 TRENDPARAM 	<- as.integer(9)
 NUGGETVAR 	<- as.integer(10)
 CRITICALPARAM 	<- as.integer(11)
 DONOTVERIFYPARAM 	<- as.integer(12)
 ONLYRETURN 	<- as.integer(13)
 FORBIDDENPARAM 	<- as.integer(14)
 UNKNOWNPARAM 	<- as.integer(15)
 VARONLYMLE 	<- as.integer(16)
 CRITONLYMLE 	<- as.integer(17)
 ONLYMLE 	<- as.integer(18)
 IGNOREPARAM 	<- as.integer(19)
 STANDARD 	<- as.integer(20)
 INCLUDENOTRETURN 	<- as.integer(21)
 INTERNALPARAMETERS 	<- as.integer(22)
 ALLPARAMETERS 	<- as.integer(23)
 NOPARAMETERS 	<- as.integer(24)


 LASTRETURNED 	<- as.integer(FORBIDDENPARAM)
 FIRSTONLYMLE 	<- as.integer(VARONLYMLE)
 LASTONLYMLE 	<- as.integer(ONLYMLE)
 LASTUSERSORTOF 	<- as.integer(FORBIDDENPARAM)
 LASTSORTOF 	<- as.integer(NOPARAMETERS)

 CircEmbed 	<- as.integer(0)
 CircEmbedCutoff 	<- as.integer(1)
 CircEmbedIntrinsic 	<- as.integer(2)
 TBM 	<- as.integer(3)
 SpectralTBM 	<- as.integer(4)
 Direct 	<- as.integer(5)
 Sequential 	<- as.integer(6)
 Trendproc 	<- as.integer(7)
 Average 	<- as.integer(8)
 Nugget 	<- as.integer(9)
 RandomCoin 	<- as.integer(10)
 Hyperplane 	<- as.integer(11)
 Specific 	<- as.integer(12)
 Nothing 	<- as.integer(13)
 Forbidden 	<- as.integer(14)


 original_model 	<- as.integer(0)
 mle_conform 	<- as.integer(1)
 all_origins 	<- as.integer(2)


 INTERNAL_PARAM 	<- "internal"
 ADD_NA 	<- "add.na"

 GETMODEL_AS_SAVED 	<- as.integer(0)
 GETMODEL_DEL_NATSC 	<- as.integer(1)
 GETMODEL_SOLVE_NATSC 	<- as.integer(2)
 GETMODEL_DEL_MLE 	<- as.integer(3)
 GETMODEL_SOLVE_MLE 	<- as.integer(4)

 MIXED_X_NAME 	<- "X"
 MIXED_BETA_NAME 	<- "beta"

 COVARIATE_C_NAME 	<- "data"
 COVARIATE_X_NAME 	<- "x"
 COVARIATE_RAW_NAME 	<- "raw"
 COVARIATE_EXTRA_DATA_NAME 	<- "extra_data"
 COVARIATE_ADDNA_NAME 	<- "addNA"
 COVARIATE_NAME_NAME 	<- "name"

 COVARIATE_DATA_NAME 	<- "data"

 CONST_A_NAME 	<- "x"

 MINMAX_PMIN 	<- as.integer(1)
 MINMAX_PMAX 	<- as.integer(2)
 MINMAX_TYPE 	<- as.integer(3)
 MINMAX_NAN 	<- as.integer(4)
 MINMAX_MIN 	<- as.integer(5)
 MINMAX_MAX 	<- as.integer(6)
 MINMAX_OMIN 	<- as.integer(7)
 MINMAX_OMAX 	<- as.integer(8)
 MINMAX_ROWS 	<- as.integer(9)
 MINMAX_COLS 	<- as.integer(10)
 MINMAX_BAYES 	<- as.integer(11)
 MINMAX_COORD 	<- as.integer(12)
 MINMAX_ENTRIES 	<- as.integer(MINMAX_COORD)

 XLIST_X 	<- as.integer(0)
 XLIST_T 	<- as.integer(1)
 XLIST_GRID 	<- as.integer(2)
 XLIST_SPATIALDIM 	<- as.integer(3)
 XLIST_TIME 	<- as.integer(4)
 XLIST_DIST 	<- as.integer(5)
 XLIST_Y 	<- as.integer(6)
 XLIST_TY 	<- as.integer(7)
 XLIST_GRIDY 	<- as.integer(8)
 XLIST_RELEVANT_ELMTS 	<- as.integer((XLIST_GRIDY+1))

 XLIST_RESTOT 	<- as.integer(9)
 XLIST_RESTOTY 	<- as.integer(10)
 XLIST_UNITS 	<- as.integer(11)
 XLIST_NEWUNITS 	<- as.integer(12)
 XLIST_RAWXIDX 	<- as.integer(13)
 XLIST_RAWSET 	<- as.integer(14)
 XLIST_ENTRIES 	<- as.integer((XLIST_NEWUNITS+1))

 PROJ_SPACE 	<- as.integer(-1)
 PROJ_TIME 	<- as.integer(-2)
 PROJECTIONS 	<- as.integer(2)

 VAR2COV_EXTREMAL 	<- as.integer(-1)
 VAR2COV_ORIGIN 	<- as.integer(-2)
 VAR2COV_CENTER 	<- as.integer(-3)
 VAR2COV_METHODS 	<- as.integer(3)

 INTERN_SHOW 	<- as.integer(2)

RC_VARIOGRAM <- VARIOGRAM 	<- as.integer(-2)
RC_MADOGRAM <- MADOGRAM 	<- as.integer(-1)
RC_COVARIANCE <- COVARIANCE 	<- as.integer(0)

RC_PSEUDOMADOGRAM <- PSEUDOMADOGRAM 	<- as.integer(1)
RC_PSEUDO <- PSEUDO 	<- as.integer(2)
 TOTAL_FCTN_TYPE 	<- as.integer(5)
 ALPHAPSEUDOMADOGRAM 	<- as.integer(3)

 EV_FFT_EV 	<- as.integer(0)
 EV_FFT_N 	<- as.integer(1)
 EV_FFT_VAR 	<- as.integer(2)
 EV_EV 	<- as.integer(0)
 EV_SDSQ 	<- as.integer(1)
 EV_N 	<- as.integer(2)

 BRP_UNIF 	<- as.integer(0)
 BRP_COV 	<- as.integer(1)
 BRP_KRIG 	<- as.integer(2)

 VAR2COV_X 	<- as.integer(0)
 VAR2COV_C 	<- as.integer(1)

 POISSON_SCATTER_OPTIM 	<- as.integer(0)
 POISSON_SCATTER_STANDARD 	<- as.integer(1)
 POISSON_SCATTER_ANY 	<- as.integer(2)
 nPOISSON_SCATTER 	<- as.integer((POISSON_SCATTER_ANY+1))
 nEQ_NAMES 	<- as.integer(6)

 NLOPTR 	<- as.integer(3)

 XSTART 	<- as.integer(0)
 XSTEP 	<- as.integer(1)
 XLENGTH 	<- as.integer(2)

 NEIGHB_MIN 	<- as.integer(0)
 NEIGHB_SPLIT 	<- as.integer(1)
 NEIGHB_MED 	<- as.integer(NEIGHB_SPLIT)
 NEIGHB_MAX 	<- as.integer(2)

 BR_NAME 	<- "brownresnick"
 EG_NAME 	<- "schlather"
 OPITZ_NAME 	<- "opitz"
 SMITH_NAME 	<- "smith"

 MAX_NA 	<- as.integer(100)



 ## from  src/AutoRandomFields.cc

RC_DOMAIN_NAMES <- DOMAIN_NAMES <-
c( "single variable","kernel","framework dependent","submodel dependent","framework or submodel dependent","parameter dependent","<keep copy>","mismatch" )


RC_OPTIMISER_NAMES <- OPTIMISER_NAMES <-
c( "optim","optimx","soma","nloptr","GenSA","minqa","pso","DEoptim" )


RC_NLOPTR_NAMES <- NLOPTR_NAMES <-
c( "NLOPT_GN_DIRECT","NLOPT_GN_DIRECT_L","NLOPT_GN_DIRECT_L_RAND","NLOPT_GN_DIRECT_NOSCAL","NLOPT_GN_DIRECT_L_NOSCAL","NLOPT_GN_DIRECT_L_RAND_NOSCAL","NLOPT_GN_ORIG_DIRECT","NLOPT_GN_ORIG_DIRECT_L","NLOPT_LN_PRAXIS","NLOPT_GN_CRS2_LM","NLOPT_LN_COBYLA","NLOPT_LN_NELDERMEAD","NLOPT_LN_SBPLX","NLOPT_LN_BOBYQA","NLOPT_GN_ISRES" )


RC_LIKELIHOOD_NAMES <- LIKELIHOOD_NAMES <-
c( "auto","full","composite","tesselation" )


RC_DUPLICATEDLOC_NAMES <- DUPLICATEDLOC_NAMES <-
c( "error","risk error","regard as repeated","scatter" )


RC_ISO_NAMES <- ISO_NAMES <-
c( "isotropic","space-isotropic","vector-isotropic","symmetric","cartesian system","gnomonic projection","orthographic projection","spherical isotropic","spherical symmetric","spherical system","earth isotropic","earth symmetric","earth system","cylinder system","non-dimension-reducing","framework dependent","submodel dependent","framework or submodel dependent","parameter dependent","<internal keep copy>","<mismatch>" )


RC_TYPE_NAMES <- TYPE_NAMES <-
c( "tail correlation","positive definite","variogram","negative definite","point-shape function","shape function","trend","distribution or shape","of manifold type","process","method for Gauss process","normed process (non-negative values with maximum value being 0 or 1)","method for Brown-Resnick process","Smith","Schlather","Poisson","PoissonGauss","distribution family","interface","mathematical operator","other type","badly defined","<same as previous>","likelihood","evaluation","mixed input","string input","<special I>","<special II>","<special III>","<special IV>" )


 NEGATIVE_NAMES <-
c( "PARAM_DEP","PREVMODEL_DEP","SUBMODEL_DEP","MISMATCH","UNSET" )


RC_MONOTONE_NAMES <- MONOTONE_NAMES <-
c( "not set","mismatch in monotonicity","submodel dependent monotonicity","previous model dependent monotonicity","parameter dependent monotonicity","not monotone","monotone","Gneiting-Schaback class","normal mixture","completely monotone","Bernstein" )


 SORT_ORIGIN_NAMES <-
c( "original model","MLE conform","all" )


 MODE_NAMES <-
c( "careless","sloppy","easygoing","normal","precise","pedantic","neurotic" )


 OUTPUTMODE_NAMES <-
c( "sp","RandomFields","geoR" )


 ANGLE_NAMES <-
c( "radians","degree" )


 REPORTCOORD_NAMES <-
c( "always","warn","important","never" )


 UNITS_NAMES <-
c( "","km","miles","<time>","<user defined>" )


 COORD_SYS_NAMES <-
c( "auto","keep","cartesian","earth","sphere","gnomonic","orthographic","coordinate system mix" )


 COORD_NAMES_GENERAL <-
c( "coords.x","coords.T" )


 CARTESIAN_SYS_NAMES <-
c( "cartesian","gnomonic","orthographic" )


 TYPEOF_PARAM_NAMES <-
c( "variance","covariance","sd","signed sd","scale","diagonal","aniso","integer","unspecified","trend","nugget","critical to estimate","never verified","internally ignored","forbidden to be estimated","unkown","variance (used internally)","critical to estimate (used internally)","only used by MLE","neither used by MLE nor returned","standard","include never returned","internal","all","none" )


 EQ_NAMES <-
c( "==","!=","<=","<",">=",">" )


 NAMES_OF_NAMES <-
c( "EQ_NAMES","ISO_NAMES","DOMAIN_NAMES","TYPE_NAMES","MONOTONE_NAMES","MODE_NAMES","OUTPUTMODE_NAMES","REPORTCOORD_NAMES","UNITS_NAMES","COORD_SYS_NAMES","CARTESIAN_SYS_NAMES","TYPEOF_PARAM_NAMES" )


 PROJECTION_NAMES <-
c( "space","time" )


 RMCOV_X <-
c( "origin","center","extremals","all" )


 FCTN_TYPE_NAMES <-
c( "Variogram","Madogram","Covariance","Pseudomadogram","Pseudovariogram" )


RC_METHOD_NAMES <- METHOD_NAMES <-
c( "circulant","cutoff","intrinsic","tbm","spectral","direct","sequential","trend","average","nugget","coins","hyperplane","specific","any method","forbidden" )


 POISSON_SCATTER_NAMES <-
c( "optimized","standard","any" )



list2RMmodel_Names <- c('iR_P', 'iR.c', 'iR.p', 'iRFcovariance', 'iRFcovmatrix', 'iRFloglikelihood', 'iRFmadogram', 'iRFpseudomadogram', 'iRFpseudovariogra', 'iRFsimulate', 'iRFvariogram', 'iRMcov', 'iRMcovariate', 'iRMfixcov', 'R.acos', 'R.acosh', 'R.asin', 'R.asinh', 'R.atan', 'R.atan2', 'R.atanh', 'R.c', 'R.cbrt', 'R.ceil', 'R.const', 'R.cos', 'R.cosh', 'R.div', 'R.erf', 'R.erfc', 'R.exp', 'R.exp2', 'R.expm1', 'R.fabs', 'R.fdim', 'R.floor', 'R.fmax', 'R.fmin', 'R.fmod', 'R.gamma', 'R.ggamma', 'R.hypot', 'R.is', 'R.lat', 'R.lgamma', 'R.log', 'R.log1p', 'R.log2', 'R.lon', 'R.minus', 'R.mult', 'R.p', 'R.plus', 'R.pow', 'R.remainder', 'R.round', 'R.sin', 'R.sinh', 'R.sqrt', 'R.tan', 'R.tanh', 'R.trunc', 'RFboxcox', 'RFcalc', 'RFcov', 'RFcovariance', 'RFcovmatrix', 'RFcrossvalidate', 'RFddistr', 'RFdistr', 'RFearth2cartesian', 'RFearth2dist', 'RFempiricalcovariance', 'RFempiricalmadogram', 'RFempiricalvariogram', 'RFfctn', 'RFfit', 'RFformula', 'RFfractaldim', 'RFgetMethodNames', 'RFgetModel', 'RFgetModelInfo', 'RFgetModelInfo_model', 'RFgetModelInfo_register', 'RFgetModelNames', 'RFgridDataFrame', 'RFgui', 'RFhessian', 'RFhurst', 'RFinterpolate', 'RFlikelihood', 'RFlinearpart', 'RFmadogram', 'RFoldstyle', 'RFpar', 'RFparameters', 'RFpdistr', 'RFplot', 'RFplotEmpVariogram', 'RFplotEmpVariogramX', 'RFplotModel', 'RFplotSimulation', 'RFplotSimulation1D', 'RFpointsDataFrame', 'RFpseudomadogram', 'RFpseudovariogram', 'RFqdistr', 'RFratiotest', 'RFrdistr', 'RFsimulate', 'RFspatialGridDataFrame', 'RFspatialPointsDataFrame', 'RFspDataFrame2conventional', 'RFspDataFrame2dataArray', 'RFvariogram', 'RM_COVARIATE', 'RM_DECLARE', 'RM_DISTR', 'RM_MATRIX', 'RM_MULT', 'RM_NUGGET', 'RM_PLUS', 'RM_TREND', 'RM_USER', 'RMangle', 'RMaskey', 'RMave', 'RMball', 'RMbcw', 'RMbernoulli', 'RMbessel', 'RMbicauchy', 'RMbigneiting', 'RMbistable', 'RMbiwm', 'RMblend', 'RMbr2bg', 'RMbr2eg', 'RMbrownresnick', 'RMbubble', 'RMcardinalsine', 'RMcauchy', 'RMcauchytbm', 'RMcauchyUnif1', 'RMcauchyUnif2', 'RMcauchyUnif3', 'RMchoquet', 'RMcircular', 'RMconstant', 'RMcov', 'RMCOV_X', 'RMcovariate', 'RMcoxisham', 'RMcubic', 'RMcurlfree', 'RMcutoff', 'RMdagum', 'RMdampedcos', 'RMdeclare', 'RMdelay', 'RMderiv', 'RMdewijsian', 'RMdivfree', 'RMeaxxa', 'RMepscauchy', 'RMetaxxa', 'RMexp', 'RMexponential', 'RMfbm', 'RMfixcov', 'RMflatpower', 'RMfractdiff', 'RMfractgauss', 'RMgauss', 'RMgaussGammalike', 'RMgaussgauss', 'RMgencauchy', 'RMgenfbm', 'RMgengneiting', 'RMgennsst', 'RMgneiting', 'RMgneitingdiff', 'RMhandcock', 'RMhyperbolic', 'RMiaco', 'RMid', 'RMidcoord', 'RMidmodel', 'RMintexp', 'RMintrinsic', 'RMkolmogorov', 'RMlatentCauchy1', 'RMlatentCauchy2', 'RMlatentCauchy3', 'RMlatentCauchy4', 'RMlgd', 'RMlsfbm', 'RMm2r', 'RMm3b', 'RMma', 'RMmastein', 'RMmatern', 'RMmatrix', 'RMmodelplus', 'RMmppplus', 'RMmps', 'RMmqam', 'RMmult', 'RMmultiquad', 'RMnatsc', 'RMnsst', 'RMnugget', 'RMparswm', 'RMparswmX', 'RMpenta', 'RMplus', 'RMpolygon', 'RMpolynome', 'RMpower', 'RMpoweredexp', 'RMprod', 'RMqam', 'RMqexp', 'RMrational', 'RMrotat', 'RMrotation', 'RMS', 'RMscale', 'RMschlather', 'RMschur', 'RMsign', 'RMsinepower', 'RMspheric', 'RMstable', 'RMstein', 'RMstp', 'RMsum', 'RMtbm', 'RMtent', 'RMtrafo', 'RMshape', 'RMshapeplus', 'RMtruncsupport', 'RMuser', 'RMvector', 'RMwave', 'RMwendland', 'RMwhittle', 'RPaverage', 'RPbernoulli', 'RPbrmixed', 'RPbrorig', 'RPbrownresnick', 'RPbrshifted', 'RPchi2', 'RPcirculant', 'RPcoins', 'RPcutoff', 'RPdirect', 'RPgauss', 'RPhyperplane', 'RPintrinsic', 'RPloggaussnormed', 'RPnugget', 'RPopitz', 'RPpoisson', 'RPschlather', 'RPsequential', 'RPsmith', 'RPspecific', 'RPspectral', 'RPt', 'RPtbm', 'RPtrend', 'RRdeterm', 'RRdistr', 'RRgauss', 'RRloc', 'RRmcmc', 'RRrectangular', 'RRspheric', 'RRunif')

list2RMmodel_oldNames <- c('#', '>', '>', '+', '*', '$', '$power', 'ave', 'shape.ave', 'bcw', 'locstatfbm', 'bessel', 'bigneiting', 'bernoulli', 'biWM', 'bistable', 'blend', 'brownresnick', 'br2bg', 'br2eg', 'bubble', 'cauchy', 'cauchyUnif1', 'cauchyUnif2', 'cauchyUnif3', 'latentCauchy1', 'latentCauchy2', 'latentCauchy3', 'latentCauchy4', 'circular', 'CDeWijsian', 'constant', 'fixcov', 'coxisham', 'cubic', 'curlfree', 'cutoff', 'dagum', 'dampedcosine', 'deriv', 'DeWijsian', 'divfree', 'epsC', 'exponential', 'Exp', 'extremalgauss', 'FD', 'flatpower', 'fractalB', 'fractgauss', 'gauss', 'gaussgauss', 'gaussGammalike', 'genB', 'gencauchy', 'bicauchy', 'gengneiting', 'gengneit_intern', 'gneiting', 'gennsst', 'gennsst_intern', 'hyperbolic', 'iacocesare', 'identity', 'kolmogorov', 'lgd1', 'mastein', 'ma1', 'ma2', 'M', 'matern', 'mqam', 'multiquad', 'natsc', 'nsst', 'nugget', 'parsWM', 'penta', 'power', 'Pow', 'prod', 'qam', 'qexponential', 'scale', 'schur', 'shift', 'sinepower', 'spherical', 'stable', 'Stein', 'steinst1', 'stp', 'shape.stp', 'tbm', 'sum', 'U', 'cov', 'vector', 'wave', 'whittle', 'missing', 'null', 'trend', 'select', 'angle', 'ball', 'covariate', 'declare', 'EAxxA', 'EtAxxA', 'idcoord', 'trafo', 'mult_inverse', 'polygon', 'rational', 'rotat', 'Rotat', 'scatter', 'sign', 'setparam', 'm2r', 'm3b', 'r3binner', 'mps', 'truncsupport', 'arcsqrt', 'determ', 'distr', 'normal', 'setDistr', 'loc', 'mcmc', 'rectangular', 'spheric', 'unif', 'MCMC_PGS', 'zhou', 'ballani', 'standardShape', '++', 'statShape', 'Simulate', 'Covariance', 'CovMatrix', 'Dummy', 'get', 'Fctn', 'Distr', 'loglikelihood', 'linearpart', 'predict', 'Pseudovariogram', 'Pseudomadogram', 'Madogram', 'Variogram', '$proc', 'plusproc', 'prodproc', 'trafoproc', 'mppplusproc', 'multproc', 'matrixproc', 'covproc', 'trend', 'average', 'coins', 'averageIntern', 'circulant', 'cutoff', 'cutoffIntern', 'intrinsic', 'intrinsIntern', 'direct', 'hyperplane', 'hyperIntern', 'nugget', 'nuggetIntern', 'sequential', 'spectral', 'spectralIntern', 'specific', 'tbm', 'tbmIntern', 'loggaussnormed', 'brorig', 'brorigIntern', 'brmixed', 'brmixedIntern', 'brshifted', 'brshiftIntern', 'brownresnick', 'binaryprocess', 'gauss.process', 'poisson', 'extremalgauss', 'extremalt', 'smith', 'chi2', 't', 'minus', 'plus', 'div', 'mult', 'const', 'p', 'c', 'is', '.asin', '.atan', '.atan2', '.cos', '.sin', '.tan', '.asinh', '.atanh', '.cosh', '.sinh', '.tanh', '.log', '.expm1', '.log1p', '.exp2', '.log2', '.hypot', '.cbrt', '.ceil', '.floor', '.fmod', '.round', '.trunc', '.erfc', '.lgamma', '.remainder', '.fdim', '.fmax', '.fmin', '.gamma', '.ggamma', '.exp', '.erf', '.fabs', '.acos', '.acosh', '.pow', '.sqrt')

rfgui1_Names <- c('RMaskey', 'RMbcw', 'RMbessel', 'RMcauchy', 'RMcauchytbm', 'RMcircular', 'RMcubic', 'RMdampedcos', 'RMepscauchy', 'RMexp', 'RMfractdiff', 'RMfractgauss', 'RMgauss', 'RMgencauchy', 'RMgengneiting', 'RMgneitingdiff', 'RMhyperbolic', 'RMlgd', 'RMparswm', 'RMpenta', 'RMqexp', 'RMspheric', 'RMstable', 'RMwave', 'RMwhittle')

rfgui2_Names <- c('RMaskey', 'RMbcw', 'RMbessel', 'RMcauchy', 'RMcircular', 'RMcubic', 'RMdampedcos', 'RMepscauchy', 'RMexp', 'RMgauss', 'RMgencauchy', 'RMgengneiting', 'RMhyperbolic', 'RMlgd', 'RMparswm', 'RMpenta', 'RMqexp', 'RMspheric', 'RMstable', 'RMwave', 'RMwhittle')
