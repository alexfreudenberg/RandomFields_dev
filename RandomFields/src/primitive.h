


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



#ifndef Primitives_H
#define Primitives_H 1

#define UNIT_EPSILON 1E-13
#include "RF.h"
#include "Coordinate_systems.h" // noetig ??

double interpolate(double y, double *stuetz, int nstuetz, int origin,
		   double lambda, int delta);



#define ANGLE_ANGLE 0
#define ANGLE_LATANGLE 1
#define ANGLE_RATIO 2
#define ANGLE_DIAG 3
void AngleMatrix(model *cov, double *A);
void kappa_Angle(int i, model *cov, int *nr, int *nc);
void Angle(double *x, int *, model *cov, double *v);
int checkAngle(model *cov);
void rangeAngle(model *cov, range_type* ra);
void invAngle(double *x, model *cov, double *v);


void Bessel(double *x, int *, model *cov, double *v);
int initBessel(model *cov, gen_storage *s);
void spectralBessel(model *cov, gen_storage *s, double *e); 
int checkBessel(model *cov);
void rangeBessel(model *cov, range_type *ra);


/* Cauchy models */
void Cauchy(double *x, int *, model *cov, double *v);
void logCauchy(double *x, int *, model *cov, double *v, double *Sign);
void TBM2Cauchy(double *x, int *, model *cov, double *v);
void DCauchy(double *x, int *, model *cov, double *v);
void DDCauchy(double *x, int *, model *cov, double *v);
void InverseCauchy(double *x, model *cov, double *v);
int checkCauchy(model *cov);
void rangeCauchy(model *cov, range_type *ra);
void coinitCauchy(model *cov, localfactstype *li);
void DrawMixCauchy(model *cov, double *random); 
double LogMixDensCauchy(double *x, double logV, model *cov);
int initCauchy(model *cov, gen_storage *s);
void doCauchy(model *cov, gen_storage *s);


/* another Cauchy model */
void Cauchytbm(double *x, int *, model *cov, double *v);
void DCauchytbm(double *x, int *, model *cov, double *v);
void rangeCauchytbm(model *cov, range_type *ra);


/* circular model */
void circular(double *x, int *, model *cov, double *v);
void Inversecircular(double *x, model *cov, double *v);
void Dcircular(double *x, int *, model *cov, double *v);
int structcircular(model *cov,  model **);


void kappaconstant(int i, model *cov, int *nr, int *nc);
void constant(double *x, int *, model *cov, double *v);
void nonstatconstant(double *x, double *y, int *, model *cov, double *v);
void rangeconstant(model *cov, int k, int i, int j, simple_range_type *range);
int checkconstant(model *cov) ;


/* coxgauss, cmp with nsst1 !! */
// see NonIsoCovFct.cc

/* cubic */
void cubic(double *x, int *, model *cov, double *v);
void Dcubic(double *x, int *, model *cov, double *v); 
void Inversecubic(double *x, model *cov, double *v);


/* cutoff */
// see Hypermodel.cc

/* dagum */
void dagum(double *x, int *, model *cov, double *v);
void Inversedagum(double *x, model *cov, double *v); 
void Ddagum(double *x, int*, model *cov, double *v);
void rangedagum(model *cov, range_type* ra);
int checkdagum(model *cov);
int initdagum(model *cov, gen_storage *stor);


/*  damped cosine -- derivative of exponential:*/
void dampedcosine(double *x, int *, model *cov, double *v);
void logdampedcosine(double *x, int *, model *cov, double *v, double *Sign);
void Inversedampedcosine(double *x, model *cov, double *v); 
void Ddampedcosine(double *x, int *, model *cov, double *v);
void rangedampedcosine(model *cov, range_type* ra);
int checkdampedcosine(model *cov);


void kappa_declare(int i, model *cov, int *nr, int *nc);
void declarefct(double *x, int *, model *cov, double *v);
void nonstat_declarefct(double *x, double *y, int *, model *cov, double *v);
void rangedeclare(model *cov, range_type *range);
int checkdeclare(model *cov);


/* De Wijsian */
void dewijsian(double *x, int *, model *cov, double *v);
void Ddewijsian(double *x, int *, model *cov, double *v);
void DDdewijsian(double *x, int *, model *cov, double *v);
void D3dewijsian(double *x, int *, model *cov, double *v);
void D4dewijsian(double *x, int *, model *cov, double *v);
void rangedewijsian(model *cov, range_type* ra);
int checkdewijsian(model *cov);
void Inversedewijsian(double *x, model *cov, double *v); 
void coinitdewijsian(model *cov, localfactstype *li);

/* De Wijsian */
void DeWijsian(double *x, int *, model *cov, double *v);
void rangeDeWijsian(model *cov, range_type* ra);
void InverseDeWijsian(double *x, model *cov, double *v); 


/* exponential model */
void exponential(double *x, int *, model *cov, double *v);
void logexponential(double *x, int *, model *cov, double *v, double *Sign);
void TBM2exponential(double *x, int *, model *cov, double *v);
void Dexponential(double *x, int *, model *cov, double *v);
void DDexponential(double *x, int *, model *cov, double *v);
void Inverseexponential(double *x, model *cov, double *v);
void inverse_log_nonstatExp(double *x, model *cov,double *left, double *right);
int initexponential(model *cov, gen_storage *s);
void spectralexponential(model *cov, gen_storage *s, double *e);
int checkexponential(model *cov);
int hyperexponential(double radius, double *center, double *rx,
		     model *cov, bool simulate, 
		     double ** Hx, double ** Hy, double ** Hr);
void coinitExp(model *cov, localfactstype *li);
void ieinitExp(model *cov, localfactstype *li);
void DrawMixExp(model *cov, double *random);
double LogMixDensExp(double *x, double logV, model *cov);
int init_exp(model *cov, gen_storage *s); 
void do_exp(model *cov, gen_storage *s);


// Brownian motion 
#define BROWN_ALPHA 0
#define BROWN_GEN_BETA 1
void fractalBrownian(double *x, int *, model *cov, double *v);
void fractalBrownian(double *x, int *, model *cov, double *v, double *Sign);
void logfractalBrownian(double *x, int *, model *cov, double *v, double *Sign);
void DfractalBrownian(double *x, int *, model *cov, double *v);
void DDfractalBrownian(double *x, int *, model *cov, double *v);
void D3fractalBrownian(double *x, int *, model *cov, double *v);
void D4fractalBrownian(double *x, int *, model *cov, double *v);
void rangefractalBrownian(model *cov, range_type* ra);
void ieinitBrownian(model *cov, localfactstype *li);
void coinitBrownian(model *cov, localfactstype *li);
void InversefractalBrownian(double *x, model *cov, double *v);
int checkfractalBrownian(model *cov);
int initfractalBrownian(model *cov, gen_storage *s); 



#define LOCALLY_BROWN_ALPHA BROWN_ALPHA
#define LOCALLY_BROWN_C (BROWN_ALPHA + 1)
int checklsfbm(model *cov);
void rangelsfbm(model VARIABLE_IS_NOT_USED *cov, range_type *range);
void lsfbm(double *x, int *, model *cov, double *v);
void Dlsfbm(double *x, int *, model *cov, double *v);
void DDlsfbm(double *x, int *, model *cov, double *v);
void D3lsfbm(double *x, int *, model *cov, double *v); 
void D4lsfbm(double *x, int *, model *cov, double *v); 
void Inverselsfbm(double *x, model *cov, double *v);
int initlsfbm(model *cov, gen_storage *s);



/* FD model */
int checkFD(model *cov);
void FD(double *x, int *, model *cov, double *v);
void rangeFD(model *cov, range_type* ra);



/* fractgauss */
void fractGauss(double *x, int *, model *cov, double *v);
void rangefractGauss(model *cov, range_type* ra);


/* Gausian model */
void Gauss(double *x, int *, model *cov, double *v);
void logGauss(double *x, int *, model *cov, double *v, double *Sign);
void DGauss(double *x, int *, model *cov, double *v);
void DDGauss(double *x, int *, model *cov, double *v);
void D3Gauss(double *x, int *, model *cov, double *v);
void D4Gauss(double *x, int *, model VARIABLE_IS_NOT_USED *cov, double *v);
void spectralGauss(model *cov, gen_storage *s, double *e);   
void DrawMixGauss(model *cov, double *random);
double LogMixDensGauss(double *x, double logV, model *cov);
void InverseGauss(double *x, model *cov, double *v);
void inverse_log_nonstatGauss(double *x, model VARIABLE_IS_NOT_USED *cov, 
			      double *left, double *right);
int struct_Gauss(model *cov, model **);
int initGauss(model *cov, gen_storage *s);
void do_Gauss(model *cov, gen_storage *s) ; 
//void getMassGauss(double *a, model *cov, double *kappas, double *m);
//void densGauss(double *x, model *cov, double *v);
//void simuGauss(model *cov, int dim, double *v);

void bcw(double *x, int *, model *cov, double *v);
//void logbcw(double *x, int *, model *cov, double *v, double *Sign);
void Dbcw(double *x, int *, model *cov, double *v);
void DDbcw(double *x, int *, model *cov, double *v);
void D3bcw(double *x, int *, model *cov, double *v);
void D4bcw(double *x, int *, model *cov, double *v);
int checkbcw(model *cov);
void rangebcw(model *cov, range_type* ra);
void coinitbcw(model *cov, localfactstype *li);
void ieinitbcw(model *cov, localfactstype *li);
void Inversebcw(double *x, model *cov, double *v);


/* generalised fractal Brownian motion */
void genBrownian(double *x, int *, model *cov, double *v);
void loggenBrownian(double *x, int *, model *cov, double *v, double *Sign);
void InversegenBrownian(double *x, model *cov, double *v) ;
int checkgenBrownian(model *cov);
void rangegenBrownian(model *cov, range_type* ra);


/* epsC */
void epsC(double *x, int *, model *cov, double *v);
void logepsC(double *x, int *, model *cov, double *v, double *Sign);
void DepsC(double *x, int *, model *cov, double *v);
void DDepsC(double *x, int *, model *cov, double *v);
int checkepsC(model *cov);
void rangeepsC(model *cov, range_type* ra);
void inverseepsC(double *x, model *cov, double *v);


/* gencauchy */
void generalisedCauchy(double *x, int *, model *cov, double *v);
void loggeneralisedCauchy(double *x, int *, model *cov, double *v, double *Sign);
void DgeneralisedCauchy(double *x, int *, model *cov, double *v);
void DDgeneralisedCauchy(double *x, int *, model *cov, double *v);
int checkgeneralisedCauchy(model *cov);
void rangegeneralisedCauchy(model *cov, range_type* ra);
void coinitgenCauchy(model *cov, localfactstype *li);
void ieinitgenCauchy(model *cov, localfactstype *li);
void InversegeneralisedCauchy(double *x, model *cov, double *v);


/* bivariate Cauchy bivariate genCauchy bivariate generalized Cauchy */
#define BICauchyalpha 0
#define BICauchybeta 1
#define BICauchyscale 2
#define BICauchyrho 3
void biCauchy (double *x, int *, model *cov, double *v);
void DbiCauchy(double *x, int *, model *cov, double *v);
void DDbiCauchy(double *x, int *, model *cov, double *v);
void D3biCauchy(double *x, int *, model *cov, double *v);
void D4biCauchy(double *x, int *, model *cov, double *v);
void kappa_biCauchy(int i, model *cov, int *nr, int *nc);
int checkbiCauchy(model *cov);
void rangebiCauchy(model VARIABLE_IS_NOT_USED *cov, range_type *range);
void coinitbiCauchy(model *cov, localfactstype *li);

/* gengneiting */
#define GENGNEITING_K 0
#define GENGNEITING_MU 1
void genGneiting(double *x, int *, model *cov, double *v);
void DgenGneiting(double *x, int *, model *cov, double *v);
void DDgenGneiting(double *x, int *, model *cov, double *v);
void rangegenGneiting(model *cov, range_type* ra);
int checkgenGneiting(model *cov);


/* Gneiting's functions -- alternative to Gaussian */
// #define Sqrt2TenD47 0.30089650263257344820 /* approcx 0.3 ?? */
void Gneiting(double *x, int *, model *cov, double *v); 
//void InverseGneiting(double *x, int *, model *cov, double *v);
void DGneiting(double *x, int *, model *cov, double *v); 
void DDGneiting(double *x, int *, model *cov, double *v); 
int checkGneiting(model *cov);
void rangeGneiting(model *cov, range_type *range);


/* hyperbolic */
void hyperbolic(double *x, int *, model *cov, double *v); 
void loghyperbolic(double *x, int *, model *cov, double *v, double *Sign); 
void Dhyperbolic(double *x, int *, model *cov, double *v); 
int checkhyperbolic(model *cov);
void rangehyperbolic(model *cov, range_type* ra);
int inithyperbolic(model *cov, gen_storage *s);



/* iaco cesare model */
void IacoCesare(double *x, int *, model *cov, double *v);
//void InverseIacoCesare(model *cov); return 1.0; } 
void rangeIacoCesare(model *cov, range_type* ra);
//int checkIacoCesare(model *cov); 

void Kolmogorov(double *x, int *, model *cov, double *v);
  int checkKolmogorov(model *cov);


/* local-global distinguisher */
void lgd1(double *x, int *, model *cov, double *v);
// void Inverselgd1(int *, model *cov); // Fehler in der Progammierung?
void Dlgd1(double *x, int *, model *cov, double *v);
void rangelgd1(model *cov, range_type* ra);
int checklgd1(model *cov);
//void Inverselgd1(model *cov, double u);

void multiquad(double *x, int *, model *cov, double *v);
void rangemultiquad(model VARIABLE_IS_NOT_USED *cov, range_type *range);

/* mastein */
// see Hypermodel.cc


/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */ 
#define WM_NU 0
#define WM_NOTINV 1
#define PARSnudiag 0
void Matern(double *x, int *, model *cov, double *v);
void NonStMatern(double *x, double *y, int *, model *cov, double *v);
void logNonStMatern(double *x, double*y, 
		    int *, model *cov, double *v, double *Sign);
void logMatern(double *x, int *, model *cov, double *v, double *Sign);
void DMatern(double *x, int *, model *cov, double *v);
void DDMatern(double *x, int *, model *cov, double *v);
void D3Matern(double *x, int *, model *cov, double *v);
void D4Matern(double *x, int *, model *cov, double *v);
int checkMatern(model *cov);
void rangeWM(model *cov, range_type* ra);
void ieinitWM(model *cov, localfactstype *li);
void coinitWM(model *cov, localfactstype *li);
double densityMatern(double *x, model *cov);
int initMatern(model *cov, gen_storage *s);
void spectralMatern(model *cov, gen_storage *s, double *e); 
void InverseMatern(double *x, model *cov, double *v);
void inversenonstatMatern(double *x, model *cov, double *left, double*right);

#define WM_LOGGAMMA 0
#define WM_GAMMA 1
double Intern_DWM(double x, double nu, double loggamma, double factor);


/* nugget effect model */
#define NUGGET_TOL 0
#define NUGGET_VDIM 1
#define NUGGET_PROC_TOL (COMMON_GAUSS + 1)
#define NUGGET_PROC_VDIM (COMMON_GAUSS + 2)
void nugget(double *x, int *, model *cov, double *v);
void nonstat_nugget(double *x, double *y, int *, model *cov, double *v);
void covmatrix_nugget(model *cov, bool, double *v);
char iscovmatrix_nugget(model *cov);
void inverse_nugget(double *x, model *cov, double *v); 
void inversenonstat_nugget(double *x, model *cov, double *left, double*right);
int check_nugget(model *cov);
void range_nugget(model *cov, range_type *range);
bool allowedDnugget(model *cov);
bool allowedInugget(model *cov);
bool setnugget(model *cov);
Types Typenugget(Types required, model *cov, isotropy_type requ_iso);


#define NULL_TYPE 0
void NullModel(double *x, int *, model *cov, double *v);
Types TypeNullModel(Types required, model *cov, isotropy_type i);
void rangeNullModel(model VARIABLE_IS_NOT_USED *cov, range_type *range);


/* penta */
void penta(double *x, int *, model *cov, double *v);
void Dpenta(double *x, int *, model *cov, double *v); 
void Inversepenta(double *x, model *cov, double *v);


/* power model */ 
void power(double *x, int *, model *cov, double *v);
void Inversepower(double *x, model *cov, double *v); 
void TBM2power(double *x, int *, model *cov, double *v);
void Dpower(double *x, int *, model *cov, double *v);
int checkpower(model *cov);
void rangepower(model *cov, range_type* ra);
int structpower(model *cov, model **newmodel);
int initpower(model *cov, gen_storage *s);


/* qexponential -- derivative of exponential */
void qexponential(double *x, int *, model *cov, double *v);
void Inverseqexponential(double *x, model *cov, double *v);
void Dqexponential(double *x, int *, model *cov, double *v);
void rangeqexponential(model *cov, range_type* ra);

void kappa_rational(int i, model *cov, int *nr, int *nc);
void rational(double *x, int *, model *cov, double *v);
int checkrational(model *cov);
void rangerational(model *cov, range_type* ra);

void sinepower(double *x, int *, model *cov, double *v);
void rangesinepower(model VARIABLE_IS_NOT_USED *cov, range_type *range);

/* spherical model */ 
void spherical(double *x, int *, model *cov, double *v);
void Inversespherical(double *x, model *cov, double *v);
void TBM2spherical(double *x, int *, model *cov, double *v);
void Dspherical(double *x, int *, model *cov, double *v);
void DDspherical(double *x, int *, model *cov, double *v);
int structspherical(model *cov, model **);
int initspherical(model *cov, gen_storage *s);
void dospherical(model *cov, gen_storage *s);


/* stable model */
void stable(double *x, int *, model *cov, double *v);
void logstable(double *x, int *, model *cov, double *v, double *Sign);
void Dstable(double *x, int *, model *cov, double *v);
void DDstable(double *x, int *, model *cov, double *v);
int checkstable(model *cov);  
void rangestable(model *cov, range_type* ra);
void coinitstable(model *cov, localfactstype *li);
void ieinitstable(model *cov, localfactstype *li);
void Inversestable(double *x, model *cov, double *v);
void inverse_log_nonstat_stable(double *x, model *cov, double *L, double *R);


/* DOUBLEISOTROPIC stable model for testing purposes only */
void stableX(double *x, int *, model *cov, double *v);
void DstableX(double *x, int *, model *cov, double *v);


// bivariate exponential or bivariate stable model
#define BIStablealpha 0
#define BIStablescale 1
#define BIStablecdiag 2
#define BIStablerho 3
#define BIStablerhored 4
#define BIStablebetared 5
#define BIStablealphadiag 6
void kappa_biStable(int i, model *cov, int *nr, int *nc);
void biStable (double *x, int *, model *cov, double *v);
void DbiStable (double *x, int *, model *cov, double *v);
void DDbiStable (double *x, int *, model *cov, double *v);
void D3biStable (double *x, int *, model *cov, double *v);
void D4biStable (double *x, int *, model *cov, double *v);
int checkbiStable(model *cov);
void rangebiStable(model *cov, range_type *range);
int initbiStable(model *cov, gen_storage *stor);
void coinitbiStable(model *cov, localfactstype *li);
sortsofparam sortof_bistable(model VARIABLE_IS_NOT_USED *cov, int k, int row,
			  int col, sort_origin origin);

/* Stein */
// see Hypermodel.cc

/* stein space-time model */
void kappaSteinST1(int i, model *cov, int *nr, int *nc);
void SteinST1(double *x, int *, model *cov, double *v);
int initSteinST1(model *cov, gen_storage *s);
void spectralSteinST1(model *cov, gen_storage *s, double *e);
void rangeSteinST1(model *cov, range_type* ra);
int checkSteinST1(model *cov);  


/* wave */
void wave(double *x, int *, model *cov, double *v);
void Inversewave(double *x, model *cov, double *v);
int initwave(model *cov, gen_storage *s);
void spectralwave(model *cov, gen_storage *s, double *e); 

/* Whittle-Matern or Whittle or Besset */ 
void Whittle(double *x, int *, model *cov, double *v);
void NonStWhittle(double *x, double *y, int *, model *cov, double *v);
void logNonStWhittle(double *x, double*y, 
		    int *, model *cov, double *v, double *Sign);
void logWhittle(double *x, int *, model *cov, double *v, double *Sign);
void TBM2Whittle(double *x, int *, model *cov, double *v);
void DWhittle(double *x, int *, model *cov, double *v);
void DDWhittle(double *x, int *, model *cov, double *v);
void D3Whittle(double *x, int *, model *cov, double *v);
void D4Whittle(double *x, int *, model *cov, double *v);
bool allowedDWM(model *cov);
bool allowedIWM(model *cov);
bool setWM(model *cov);
Types TypeWM(Types required, model *cov, isotropy_type required_iso);
int initWM(model *cov, gen_storage *stor);
int checkWM(model *cov);
double densityWhittle(double *x, model *cov);
int initWhittle(model *cov, gen_storage *s);
void spectralWhittle(model *cov, gen_storage *s, double *e); 
void DrawMixWM(model *cov, double *random);
double LogMixDensW(double *x, double logV, model *cov);
void InverseWhittle(double *x, model *cov, double *v);
void inversenonstatWhittle(double *x, model *cov, double *left, double*right);

double WM(double x, double nu, double factor);
double logWM(double x, double nu1, double nu2, double factor);
double DWM(double x, double nu, double factor);

#define GNEITING_K GENGNEITING_K    // important to keep !
#define GNEITING_MU 1
#define GNEITING_S 2
#define GNEITING_SRED 3
#define GNEITING_GAMMA 4
#define GNEITING_CDIAG 5
#define GNEITING_RHORED 6
#define GNEITING_C 7
void kappa_biGneiting(int i, model *cov, int *nr, int *nc);
void biGneiting(double *x, int *, model *cov, double *v);
void DbiGneiting(double *x, int *, model *cov, double *v);
void DDbiGneiting(double *x, int *, model *cov, double *v);
int checkbiGneiting(model *cov);
void rangebiGneiting(model *cov, range_type* ra);
int initbiGneiting(model *cov, gen_storage *s);
sortsofparam sortof_biGneiting(model VARIABLE_IS_NOT_USED *cov, int k, int row,
			       int col, sort_origin origin);


/* User defined model */
void kappaUser(int i, model *cov, int *nr, int *nc);
void User(double *x, int *, model *cov, double *v);
void nonstatUser(double *x, double *y, int *, model *cov, double *v);
void DUser(double *x, int *, model *cov, double *v);
void DDUser(double *x, int *, model *cov, double *v);
int checkUser(model *cov);
void rangeUser(model *cov, range_type *ra);
Types TypeUser(Types required, model *cov, isotropy_type i);
bool allowedDuser(model *cov);
bool allowedIuser(model *cov);
bool setUser(model *cov);


#define BInudiag 0
#define BInured 1
#define BInu 2
#define BIs 3
#define BIcdiag 4
#define BIrhored 5
#define BIc 6
#define BInotinvnu 7
void kappa_biWM(int i, model *cov, int *nr, int *nc);
void biWM2(double *x, int *, model *cov, double *v);
void biWM2D(double *x, int *, model *cov, double *v);
void biWM2DD(double *x, int *, model *cov, double *v);
void biWM2D3(double *x, int *, model *cov, double *v);
void biWM2D4(double *x, int *, model *cov, double *v);
int checkbiWM2(model *cov);
void rangebiWM2(model *cov, range_type* ra);
int initbiWM2(model *cov, gen_storage *s);
void coinitbiWM2(model *cov, localfactstype *li);
sortsofparam sortof_biwm2(model VARIABLE_IS_NOT_USED *cov, int k, int row,
			  int col, sort_origin origin);
void rangebiWM2(model *cov, range_type* ra);


void kappa_parsWM(int i, model *cov, int *nr, int *nc);
void parsWM(double *x, int *, model *cov, double *v);
void parsWMD(double *x, int *, model *cov, double *v);
int checkparsWM(model *cov);
void rangeparsWM(model *cov, range_type* ra);
int initparsWM(model *cov, gen_storage *s);

extern double WhittleUpperNu[Nothing + 1];
#endif /* Primitives_H*/
