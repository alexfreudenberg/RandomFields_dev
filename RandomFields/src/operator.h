


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


#ifndef Operators_H
#define Operators_H 1

#include "primitive.h"
#include "families.h"
#include "Coordinate_systems.h"

// ----------------------------------------------------------------------
// gneiting.cc, operator.cc, plusmalS.cc


void nonstatscale(double *x, double *y, int *, model *cov, double *v);  
int checkscale(model *cov);

#define BUBBLE_COV 0
#define BUBBLE_SCALE 1
#define BUBBLE_Z 0
#define BUBBLE_WEIGHT 1
#define BUBBLE_MINSCALE 2 // must be decreasing-- will be turned into tau=1/s^2!
#define BUBBLE_BARYCENTRE 3
void kappabubble(int i, model *cov, int *nr, int *nc);
void nonstatbubble(double *x, double *y, int *, model *cov, double *v);  
int checkbubble(model *cov);
void rangebubble(model VARIABLE_IS_NOT_USED *cov, range_type *range);

#define BLEND_MULTI 0
#define BLEND_BLEND 1
#define BLEND_THRES 0
void kappablend(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc);
void nonstatblend(double *x, double *y, int *, model *cov, double *v);  
int checkblend(model *cov);
void rangeblend(model VARIABLE_IS_NOT_USED *cov, range_type *range);


#define AVERAGE_YPHASE 0
#define AVERAGE_YFREQ 1
void kappa_ave(int i, model *cov, int *nr, int *nc);
void ave(double *x, int *, model *cov, double *v);
int checkave(model *cov);
void rangeave(model *cov, range_type* ra);
//void sd_standard(mpp_storage *s, model *cov);
int structAve(model *cov, model **newmodel);

int check_shapeave(model *cov);
int init_shapeave(model *cov, gen_storage *s);
void do_shapeave(model *cov, gen_storage *s);
void logshapeave(double *x, int *, model *cov, double *v, double *Sign);

void kappavariogram2cov(int i, model *cov, int *nr, int *nc);
void variogram2cov(double *x, double *y, int *, model *cov, double *v);
int checkvariogram2cov(model *cov);
void rangevariogram2cov(model VARIABLE_IS_NOT_USED *cov, range_type *range);


#define PLUS_TREND 0
void plus(double *x, int *, model *cov, double *v);
void nonstat_plus(double *x,  double *y, int *, model *cov, double *v);
void covmatrix_plus(model *cov, bool,double *v);
char iscovmatrix_plus(model *cov);
int checkplus(model *cov);
bool allowedDplus(model *cov);
bool allowedIplus(model *cov);
Types Typeplus(Types required, model *cov, isotropy_type i);
void Dplus(double *x, int *, model *cov, double *v);
void DDplus(double *x, int *, model *cov, double *v);
void spectralplus(model *cov, gen_storage *s, double *e);
int structplus(model *cov, model **newmodel);
int initplus(model *cov, gen_storage *s);
void doplus(model *cov, gen_storage *s);
void rangeplus(model *cov, range_type *range);


void mal(double *x, int *, model *cov, double *v);
void nonstat_mal(double *x,  double *y, int *, model *cov, double *v);
void logmal(double *x, int *, model *cov, double *v, double *Sign);
void nonstat_logmal(double *x, double *y, int *, model *cov, double *v, 
		   double *Sign);
int checkmal(model *cov);
Types Typemal(Types required, model *cov, isotropy_type i);
void Dmal(double *x, int *, model *cov, double *v);
int structmal(model *cov, model VARIABLE_IS_NOT_USED **newmodel);
int initmal(model *cov, gen_storage *s);
void domal(model *cov, gen_storage *s);

#define M_M 0
#define M_VDIM 1
void kappaM(int i, model *cov, int *nr, int *nc);
void Matrix(double *x, int *, model *cov, double *v);
void nonstatM(double *x, double *y, int *, model *cov, double *v);
int checkM(model *cov);
bool allowedDM(model *cov);
bool allowedIM(model *cov);
void rangeM(model *cov, int k, int i, int j, simple_range_type *range);
sortsofparam sortof_M(model *cov, int k, int row, int col,
		      sort_origin original);
Types TypeM(Types required, model *cov, isotropy_type i);
int initM(model *cov, gen_storage *s);
bool setM(model *cov);


void kappaSchur(int i, model *cov, int *nr, int *nc);
void Schur(double *x, int *, model *cov, double *v);
void DSchur(double *x, int *, model *cov, double *v);
void D2Schur(double *x, int *, model *cov, double *v);
void D3Schur(double *x, int *, model *cov, double *v);
void D4Schur(double *x, int *, model *cov, double *v);
void nonstatSchur(double *x, double *y, int *, model *cov, double *v);
int checkSchur(model *cov);
void rangeSchur(model *cov, int k, int i, int j,
		      simple_range_type *range);

void Id(double *x, int *, model *cov, double *v);
void nonstatId(double *x, double *y, int *, model *cov, double *v);
void logId(double *x, int *, model *cov, double *v, double *Sign);
void inverseId(double *x, model *cov, double *v);
void inversenonstatId(double *x, model *cov, double *left, double*right);
void inverse_log_nonstatId(double *v, model *cov, double *x, double *y);
int checkId(model *cov);
void DId(double *x, int *, model *cov, double *v);
void DDId(double *x, int *, model *cov, double *v);
void TBM2Id(double *x, int *, model *cov, double *v);
int initId(model *cov, gen_storage *s);
void spectralId(model *cov, gen_storage *s, double *e);
void coinitId(model *cov, localfactstype *li);
void ieinitId(model *cov, localfactstype *li);
void rangeId(model *cov, range_type* ra); 
Types TypeId(Types required, model *cov, isotropy_type required_iso);


#define DOLLAR_SUB 0
#define DVAR 0
#define DSCALE 1
#define DANISO 2 /* internal */
#define DAUSER 3 /* user defined */
#define DPROJ 4
#define DMAX DPROJ
#define PPROJ cov->Sdollar->proj
#define Nproj cov->Sdollar->nproj

void kappaS(int i, model *cov, int *nr, int *nc);
void Siso(double *x, int *, model *cov, double *v);
void logSiso(double *x, int *, model *cov, double *v, double * Sign);
void Sstat(double *x, int *, model *cov, double *v);
void covmatrixS(model *cov, bool, double *v);
char iscovmatrixS(model *cov);
void logSstat(double *x, int *, model *cov, double *v, double * Sign);
void DS(double *x, int *, model *cov, double *v);
void DDS(double *x, int *, model *cov, double *v);
void D3S(double *x, int *, model *cov, double *v);
void D4S(double *x, int *, model *cov, double *v);
void tbm2S(double *x, int *, model *cov, double *v);
void nonstatS(double *x, double *y, int *, model *cov, double *v);  
void nonstat_logS(double *x, double *y, int *, model *cov, double *v, double *);  
int checkS(model *cov);
bool allowedDS(model *cov); // returns true when everything is OK, eg RMconst
bool allowedIS(model *cov); // dito
Types TypeS(Types required, model *cov, isotropy_type i);
void rangeS(model *cov, range_type* ra);
void spectralS(model *cov, gen_storage *s, double *e);
void coinitS(model *cov, localfactstype *li);
void ieinitS(model *cov, localfactstype *li);
void inverseS(double *x, model *cov, double *v);
void inversenonstatS(double *x, model *cov, double *left, double*right);
void inverse_log_nonstatS(double *v, model *cov, double *x, double *y);
void nablaS(double *x, int *, model *cov, double *v);
void hessS(double *x, int *, model *cov, double *v);
int structS(model *cov, model **newmodel);
int initS(model *cov, gen_storage *s);
void doS(model *cov, gen_storage *s);
void ScaleDollarToLoc(model *to, model *from, int depth);
bool ScaleOnly(model *cov);
bool isDollar(model *cov);





void binary(double *x, int *, model *cov, double *v);
int checkbinary(model *cov);
void rangebinary(model *cov, range_type *range);

void brownresnick(double *x, int *, model *cov, double *v);
int checkbrownresnick(model *cov);
void Dbrownresnick(double *x, int *, model *cov, double *v);
void DDbrownresnick(double *x, int *, model *cov, double *v);
void D3brownresnick(double *x, int *, model *cov, double *v);
int struct_brownresnick(model *cov, model **newmodel);
int init_brownresnick(model *cov, gen_storage *s);
void do_brownresnick(model *cov, gen_storage *s);

void BR2BG(double *x, int *, model *cov, double *v);
int check_BR2BG(model *cov);

void BR2EG(double *x, int *, model *cov, double *v);
int check_BR2EG(model *cov);


void kappa_cox(int i, model *cov, int *nr, int *nc);
void cox(double *x, int *, model *cov, double *v);
void rangecox(model *cov, int k, int i, int j, simple_range_type *range);
int checkcox(model *cov);
int initcox(model *cov, gen_storage *s);
void spectralcox(model *cov, gen_storage *s, double *e); 
void coxnabla(double *x, int *, model *cov, double *v);
void coxhess(double *x, int *, model *cov, double *v);

void extremalgaussian(double *x, int *, model *cov, double *v);
int check_extremalgaussian(model *cov);

void binaryGauss(double *x, int *, model *cov, double *v);
int check_binaryGauss(model *cov);

void extrgauss(double *x, int *, model *cov, double *v);
  int check_extrgauss(model *cov);


void NonStWM(double *x, double *y, int *, model *cov, double *v);
int checkNonStWM(model *cov); 
void rangeNonStWM(model *cov, range_type* ra); 
//void DrawMixNonStWM(model *cov, double *random);
//double LogMixWeightNonStWM(double *x, double logV, model *cov);


void lp(double *x, int *, model *cov, double *v);
int checklp(model *cov);
void rangelp(model *cov, range_type* ra);

void MaStein(double *x, int *, model *cov, double *v);
int check_MaStein(model *cov);
void range_MaStein(model *cov, range_type* ra);

/* nsst */
/* Tilmann Gneiting's space time models, part I */
#define NSST_DELTA 0
void nsst(double *x, int *, model *cov, double *v);
void Dnsst(double *x, int *, model *cov, double *v);
void TBM2nsst(double *x, int *, model *cov, double *v);
int checknsst(model *cov);
void rangensst(model *cov, range_type* ra);

void gennsst(double *x, int *, model *cov, double *v);
void nonstatgennsst(double *x, double *y, int *, model *cov, double *v);
int checkgennsst(model *cov);
bool allowedDgennsst(model *cov);
bool allowedIgennsst(model *cov);
void rangegennsst(model VARIABLE_IS_NOT_USED *cov, range_type* ra);

void kappa_gennsst_intern(int i, model *cov, int *nr, int *nc);
void gennsst_intern(double *x, int *, model *cov, double *v);
int checkgennsst_intern(model *cov);
void range_gennsst_intern(model *cov, range_type* range);




#define STP_S 0
#define STP_Z 1
#define STP_M 2
void kappa_stp(int i, model *cov, int *nr, int *nc);
void stp(double *x,  double *y, int *, model *cov, double *v);
int checkstp(model *cov);
void rangestp(model *cov, range_type* ra);
int structStp(model *cov, model **newmodel);


int check_shapestp(model *cov);
int init_shapestp(model *cov, gen_storage *s);
void do_shapestp(model *cov, gen_storage *s);
void logshapestp(double *x, double *u, int *, model *cov, double *v, double *);


/* undefined model -- for technical reasons only */
//void rangeundefined(model *cov, range_type* ra);
//int checkundefined(model *cov);
int checkOK(model *cov);

void kappashift(int i, model *cov, int *nr, int *nc);
void shift(double *x, int *, model *cov, double *v);
int checkshift(model *cov);
void rangeshift(model *cov, range_type* ra);

void vector(double *x, int *, model *cov, double *v);
void vectorAniso(double *x, int *, model *cov, double *v);
int checkvector(model *cov);
void rangevector(model *cov, range_type* ra);

void kappadivcurl(int i, model *cov, int *nr, int *nc);
void diverge(double *x, int *, model *cov, double *v);
void curl(double *x, int *, model *cov, double *v);
int checkdivcurl(model *cov);
void rangedivcurl(model *cov, range_type* ra);
//void rangedivcurl(model *cov, range_type* ra);


void Exp(double *x, int *, model *cov, double *v);
void DExp(double *x, int *, model *cov, double *v);
void DDExp(double *x, int *, model *cov, double *v);
int checkExp(model *cov);
void rangeExp(model *cov, range_type *range);

void ma1(double *x, int *, model *cov, double *v);
int checkma1(model *cov);
void rangema1(model *cov, range_type* ra);
void ma2(double *x, int *, model *cov, double *v);
int checkma2(model *cov);


void natsc(double *x, int *, model *cov, double *v);
void Dnatsc(double *x, int *, model *cov, double *v);
void DDnatsc(double *x, int *, model *cov, double *v);
void Inversenatsc(double *x, model *cov, double *v);
void coinitnatsc(model *cov, localfactstype *li);
void ieinitnatsc(model *cov, localfactstype *li);
void tbm2natsc(double *x, int *, model *cov, double *v);
int checknatsc(model *cov);
int initnatsc(model *cov);
void spectralnatsc(model *cov, gen_storage *s, double *e);
int initnatsc(model *cov, gen_storage *s);
void donatsc(model *cov, gen_storage *s);

#define POW_ALPHA 0
void Pow(double *x, int *, model *cov, double *v);
void DPow(double *x, int *, model *cov, double *v);
void DDPow(double *x, int *, model *cov, double *v);
int checkPow(model *cov);
void rangePow(model *cov, range_type* ra); 
void InversePow(double *x, model *cov, double *v);
void inversenonstatPow(double *x, model *cov, double *left, double*right);
int initPow(model *cov, gen_storage *s);



#define POWVAR 0
#define POWSCALE 1
#define POWPOWER 2
#define POW_SUB 0
void PowS(double *x, int *, model *cov, double *v);
void logPowS(double *x, int *, model *cov, double *v, double *Sign);
void nonstatPowS(double *x, double *y, int *, model *cov, double *v);
void nonstat_logPowS(double *x, double *y, int *, model *cov, double *v, 
		 double *Sign);
void inversePowS(double *x, model *cov, double *v) ;
void inversenonstatPowS(double *x, model *cov, double *left, double *right);
int checkPowS(model *cov) ;
Types TypePowS(Types required, model *cov, isotropy_type i) ;
void rangePowS(model *cov, range_type* range);
void PowScaleToLoc(model *to, model *from, int depth);
int structPowS(model *cov, model **newmodel) ;
int initPowS(model *cov, gen_storage *s);
void doPowS(model *cov, gen_storage *s);



#define QAM_THETA 0
void kappaqam(int i, model *cov, int *nr, int *nc);
void qam(double *x, int *, model *cov, double *v);
int checkqam(model *cov);
void rangeqam(model *cov, range_type* ra);

void kappamqam(int i, model *cov, int *nr, int *nc);
void mqam(double *x, int *, model *cov, double *v);
int checkmqam(model *cov);
void rangemqam(model *cov, range_type* ra);


#define SELECT_SUBNR 0
void select(double *x, int *, model *cov, double *v);
void covmatrix_select(model *cov, bool, double *v);
char iscovmatrix_select(model *cov);
int checkselect(model *cov);
bool allowedDselect(model *covD);
bool allowedIselect(model *cov);
void rangeselect(model *cov, range_type *range);

#define TBMOP_FULLDIM 0
#define TBMOP_TBMDIM 1
#define TBMOP_LAYERS 2
#define TBMOP_COV 0
void tbm(double *x, int *, model *cov, double *v);
void Dtbm(double *x, int *, model *cov, double *v);
int checktbmop(model *cov);
void rangetbm_common(model *cov, range_type *range, bool tbmop);
void rangetbmop(model *cov, range_type* ra);
Types Typetbm(Types required, model *cov, isotropy_type i);
bool allowedItbm(model *cov);
bool settbm(model *cov);



#define pLOC_DIAM 0 // parameter p
#define pLOC_A 1 // always last of pLOC_*
#define pLOC_R pLOC_A // always last of pLOC_*


#define LOCAL_R 0 
#define LOCAL_CUTOFF_MSG (LOCAL_R + 1)

#define CUTOFF_B (LOCAL_CUTOFF_MSG + 1)
#define CUTOFF_ASQRTR (CUTOFF_B + 1)
#define CUTOFF_THEOR (CUTOFF_ASQRTR + 1) // muss immer == 4 sein
#define CUTOFF_CONSTANT (CUTOFF_THEOR + 1)


#define CUTOFF_MAX (CUTOFF_CONSTANT + 1)   /* size of vector q */


#define CUTOFF_R (CUTOFF_MAX + 1)
#define CUTOFF_CUBE_A (CUTOFF_R + 1)
#define CUTOFF_CUBE_B (CUTOFF_CUBE_A + 1)
#define CUTOFF_CUBE_C (CUTOFF_CUBE_B + 1)
//#define CUTOFF_CONSTANT (CUTOFF_CUBE_C + 1)

#define CUTOFF_CUBE_N (CUTOFF_CUBE_C + 1)
#define CUTOFF_CUBE_M (CUTOFF_CUBE_N + 1)
#define CUTOFF_CUBE_L (CUTOFF_CUBE_M + 1)

#define CUTOFF_MULTIVARIATE_MAX (CUTOFF_CUBE_L + 1)


#define CUTOFF_THIRD_CONDITION 3
void co(double *x, int *, model *cov, double *v);
int check_co(model *cov);
bool alternativeparam_co(model *cov);
void range_co(model *cov, range_type* ra);


#define INTRINSIC_A0 (LOCAL_CUTOFF_MSG + 1)
#define INTRINSIC_A2 (INTRINSIC_A0 + 1)
#define INTRINSIC_B (INTRINSIC_A2 + 1)
#define INTRINSIC_MAX (INTRINSIC_B + 1) /* size of vector q */
#define LOCAL_MAX (INTRINSIC_MAX + CUTOFF_MAX) /* sicher ist sicher */
void Stein(double *x, int *, model *cov, double *v);
int check_Stein(model *cov);
bool alternativeparam_Stein(model *cov);
void range_Stein(model *cov, range_type* ra);

void strokorb(double *x, int *, model *cov, double *v);
int checkstrokorb(model *cov);
int init_strokorb(model *cov, gen_storage *s);
void do_strokorb(model *cov, gen_storage *s);


//void strokorbBall(double *x, model *cov, double *v);
int checkstrokorbBall(model *cov);
int struct_strokorbBall(model *cov, model **newmodel);  
//int init_strokorbBall(model *cov, gen_storage *s);
//void do_strokorbBall(model *cov, gen_storage *s);
void rangestrokorbball(model  VARIABLE_IS_NOT_USED *cov, range_type *range);


void strokorbBallInner(double *x, int *, model *cov, double *v);
int check_strokorbBallInner(model *cov);
void range_strokorbBallInner(model *cov, range_type *range);
int init_strokorbBallInner(model *cov, gen_storage *s);
void do_strokorbBallInner(model *cov, gen_storage *s);


void strokorbPoly(double *x, int *, model *cov, double *v);
int checkstrokorbPoly(model *cov);
int struct_strokorbPoly(model *cov, model **newmodel); 


void mult_inverse(double *x, int *, model *cov, double *v);
void nonstat_mult_inverse(double *x, double *y, int *, model *cov, double *v);
int checkmult_inverse(model *cov);


void addSetParam(model **newmodel, model * remote, 
		 param_set_fct set,  bool performdo, int variant);
void addSetDistr(model **newmodel, model * remote, 
		 param_set_fct set,  bool performdo, int variant);
void setparamStat(double *x, int *, model *cov, double *v);
void nonstat_setparam(double *x,  double *y, int *, model *cov, double *v);
void covmatrix_setparam(model *cov, bool, double *v);
char iscovmatrix_setparam(model *cov);
int checksetparam(model *cov);
void range_setparam(model VARIABLE_IS_NOT_USED *cov, range_type *range);
void Inverse_setparam(double *v, model *cov, double *x);
void inversenonstat_setparam(double *v, model *cov, double *x, double *y);
Types Typesetparam(Types required, model *cov, isotropy_type i);
void Dsetparam(double *x, int *, model *cov, double *v);
void DDsetparam(double *x, int *, model *cov, double *v);
void D3setparam(double *x, int *, model *cov, double *v);
void D4setparam(double *x, int *, model *cov, double *v);
void spectralsetparam(model *cov, gen_storage *s, double *e);
int initsetparam(model *cov, gen_storage *s);
void dosetparam(model *cov, gen_storage *s);

#define SCATTER_STEP 0
#define SCATTER_MAX 1
void kappaScatter(int i, model *cov, int *nr, int *nc);
void Scatter(double *xx, int *, model *cov, double *v);
int checkScatter(model *cov);
void rangeScatter(model *cov, range_type *range); 
int struct_scatter(model *cov, model **newmodel);
int init_scatter(model *cov, gen_storage *s);
void do_scatter(model *cov, gen_storage *s);


void oesting(double *x, int *, model *cov, double *v);
void Doesting(double *x, int *, model *cov, double *v);
void DDoesting(double *x, int *, model *cov, double *v); 
int checkoesting(model *cov);
int initoesting(model *cov, gen_storage VARIABLE_IS_NOT_USED *s);
void rangeoesting(model *cov, range_type *range);

void idcoord(double *x, int *, model *cov, double *v);
int checkidcoord(model *cov);

#define TRAFO_ISO 0
void kappatrafo(int i, model *cov, int *nr, int *nc);
void trafo(double *x, int *, model *cov, double *v);  
void logtrafo(double *x, int *, model *cov, double *v, 
		     double *Sign);  
void nonstattrafo(double *x, double *y, int *, model *cov, double *v);  
void nonstat_logtrafo(double *x, double *y, int *, model *cov, double *v, 
		     double *Sign);
bool allowedDtrafo(model *cov);
bool allowedItrafo(model *cov);
bool settrafo(model *cov);
int checktrafo(model *cov);
void rangetrafo(model VARIABLE_IS_NOT_USED *cov, range_type *range);
Types Typetrafo(Types required, model *cov, isotropy_type i);
#define PRODPROC_RANDOM 0

void nonstat_prod(double *x, double *y, int *, model *cov, double *v);
int checkprod(model *cov);

void nonstatsum(double *x, double *y, int *, model *cov, double *v);
int checksum(model *cov);

void kappaderivative(int i, model *cov, int *nr, int *nc);
void derivative(double *x, int *, model *cov, double *w);
int checkderivative(model *cov);
void rangederivative(model *cov, range_type *range);


#endif /* Operators_H */
 
