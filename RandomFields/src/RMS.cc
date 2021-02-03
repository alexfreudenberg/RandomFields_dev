/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of auxiliary correlation functions 

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is gno error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nonsta     tionary models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc

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
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


//#include <Rmath.h>
//#include <stdio.h> 
//#include <R_ext/Lapack.h>
#include "questions.h"
#include "operator.h"
#include "Processes.h"
#include "startGetNset.h"
#include "variogramAndCo.h"
#include "rf_interfaces.h"
#include "shape.h"

void kappaS(int i, model *cov, int *nr, int *nc){
  switch(i) {
  case DVAR : case DSCALE :
    *nr = *nc = 1;
    break;
  case DANISO :
    *nr = OWNTOTALXDIM;
    *nc = SIZE_NOT_DETERMINED;
    break;
  case DAUSER :
    *nr = SIZE_NOT_DETERMINED;
    *nc = OWNTOTALXDIM;
    break;
  case DPROJ : 
    *nr = SIZE_NOT_DETERMINED;
    *nc = 1;
    break;
  default : *nr = *nc = OUT_OF_RANGE; 
  }
}

// simple transformations as variance, scale, anisotropy matrix, etc.  
void Siso(double *x, int *info, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  double y,
    *aniso=P(DANISO),
    *scale =P(DSCALE),
    var = P0(DVAR);
  assert(cov->Sdollar->simplevar);
  
  y = *x;
  if (aniso != NULL) y = FABS(y * aniso[0]);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] : (y==0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  // letzteres darf nur passieren wenn dim = 1!!
  COV(&y, info, next, v);

  for (i=0; i<vdimSq; v[i++] *= var); 
}
  

// simple transformations as variance, scale, anisotropy matrix, etc.  
void logSiso(double *x, int *info, model *cov, double *v, double *Sign){
  model *next = cov->sub[DOLLAR_SUB];
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  double y, 
    *aniso=P(DANISO),
    *scale =P(DSCALE),
    logvar = LOG(P0(DVAR));
   assert(cov->Sdollar->simplevar);

  y = *x;
  if (aniso != NULL) y = FABS(y * aniso[0]);

  if (scale != NULL) 
    y = scale[0]>0.0 ? y / scale[0] 
      : (y == 0.0 && scale[0]==0.0) ? 0.0 : RF_INF;
      
  LOGCOV(&y, info, next, v, Sign);
  for (i=0; i<vdimSq; v[i++] += logvar); 
  
 
}
 
void Sstat(double *x, int *info, model *cov, double *v){
  logSstat(x, info, cov, v, NULL);
}




void logSstat(double *x, int *info, model *cov, double *v, double *Sign){
  assert(cov->kappasub[DSCALE] == NULL && 
	 (cov->kappasub[DANISO]==NULL || 
	  DefList[MODELNR(cov->kappasub[DANISO])].check==checkAngle)&& 
	 (cov->kappasub[DAUSER]==NULL || 
	  DefList[MODELNR(cov->kappasub[DAUSER])].check==checkAngle));
  model *next = cov->sub[DOLLAR_SUB];
    //  *Aniso = cov->kappasub[DAUSER],
    // *Scale = cov->kappasub[DSCALE];
  double 
    *scale =P(DSCALE), 
    *aniso=P(DANISO);
  int i,
    nproj = Nproj,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  bool notdelete = false;
  TALLOC_DOUBLE(z);

  if (nproj > 0) {
    int *proj = PPROJ; 
    TALLOC_GLOBAL_X1(z, nproj);
    if (scale == NULL || scale[0] > 0.0) {
      if (scale == NULL)  for (i=0; i<nproj; i++) z[i] = x[proj[i]];
      else {
	double invscale = 1.0 / scale[0];
	for (i=0; i<nproj; i++) {
	  z[i] = invscale * x[proj[i]];
	}
      }
    } else {
      // projection and aniso may not be given at the same time
      for (i=0; i<nproj; i++)
	z[i] = (x[proj[i]] == 0 && scale[0] == 0) ? 0.0 : RF_INF;
    } 
    //  } else if (Aniso != NULL) {
    //    int dim = Aniso->vdim[0];
    //    A LLOC_DOLLAR(z, dim);
    //    FCTN(x, int *info, Aniso, z);
    //    z1 = z;
    //  } else if (Scale != NULL) {
    //    int dim = Aniso->vdim[0];
    //    A LLOC_DOLLAR(z, dim);
    //    FCTN(x, int *info, Aniso, z);
    //    z1 = z;
  } else if ((notdelete = aniso==NULL && (scale == NULL || scale[0] == 1.0))) {
    z = x;
  } else {
    int xdimown = OWNTOTALXDIM;
    double *xz;
    TALLOC_GLOBAL_X1(z, xdimown);
    if (aniso!=NULL) {
      xA(x, aniso, cov->nrow[DANISO], cov->ncol[DANISO], z);
      xz = z;
    } else xz = x;    
    if (scale != NULL) {
      double s = scale[0];
      if (s > 0.0) {
	double invscale = 1.0 / s;
	for (i=0; i < xdimown; i++) z[i] = invscale * xz[i];
      } else {
	for (i=0; i < xdimown; i++)
	  z[i] = (xz[i] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
      }
    }
  }

  double var;
  if (cov->Sdollar->simplevar) {
    var = P0(DVAR);
  } else {
    FCTN(x, info, cov->kappasub[DVAR], &var);
  }


  if (Sign==NULL) {
    COV(z, info, next, v);
    for (i=0; i<vdimSq; v[i++] *= var); 
  } else {
    LOGCOV(z, info, next, v, Sign);
    var = LOG(var);
    for (i=0; i<vdimSq; v[i++] += var); 
  }

  if (!notdelete) {
    FREE_TALLOC(z);
  }

}

void nonstatS(double *x, double *y, int *info, model *cov, double *v){
  nonstat_logS(x, y, info, cov, v, NULL);
}

void nonstat_logS(double *x, double *y, int *info,
		  model *cov, double *v, double *Sign){
  model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE];
  double 
    s1 = RF_NA, s2 = RF_NA, smeanSq=RF_NA,
    *scale =P(DSCALE),
    *aniso=P(DANISO);
  int i,
    *infoY = info + INFOS_PER_COORD,
    xdimown = OWNTOTALXDIM,   
    dimz = xdimown,
    //   nr = NA_INTEGER,
    nproj = Nproj,
    vdim = VDIM0,
    vdimSq = vdim * vdim;
  bool notdelete = false;
  TALLOC_DOUBLE(z1);
  TALLOC_DOUBLE(z2);

  //  printf("%ld %ld\n", z1, z2);
  
  if (nproj > 0) {
    int *proj = PPROJ;
    TALLOC_GLOBAL_X1(z1, nproj);
    TALLOC_GLOBAL_X2(z2, nproj);
    if (scale==NULL || scale[0] > 0.0) {
      double invscale = scale==NULL ? 1.0 :  1.0 / scale[0];
      for (i=0; i<nproj; i++) {
	z1[i] = invscale * x[proj[i]];
	z2[i] = invscale * y[proj[i]];	
      }
    } else {
      double s = scale[0]; // kann auch negativ sein ...
      for (i=0; i<nproj; i++) {
	z1[i] = (x[proj[i]] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
	z2[i] = (y[proj[i]] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
      }
    }
  } else if (Aniso != NULL) {
    dimz =  Aniso->vdim[0];
    TALLOC_GLOBAL_X1(z1, dimz);
    TALLOC_GLOBAL_X2(z2, dimz);
    FCTN(x, info, Aniso, z1);
    FCTN(y, infoY, Aniso, z2);
  } else if (Scale != NULL && !isnowRandom(Scale)) {// see also code below variance
    double s;    
    TALLOC_GLOBAL_X1(z1, dimz);
    TALLOC_GLOBAL_X2(z2, dimz);     
    FCTN(x, info, Scale, &s1);
    FCTN(y, infoY, Scale, &s2);
    
    if (s1 <= 0.0 || s2  <= 0.0)
      ERR1("'%.50s' must be a positive function", KNAME(DSCALE));
    //  if (x[0]!=1.0) printf("\n%.50s x=%10g,%10g %10g y=%10g,%10g %10g; %d\n", NAME(Scale), x[0], x[1], s1, y[0], y[1], s2, xdimown);
    smeanSq = 0.5 * (s1 * s1 + s2 * s2);
    //  printf("s_x[0]=%10g %10g s=%10g %10g %10g ", x[0], x[1], s1, s2, smeanSq);
    s = SQRT(smeanSq);
    for (i=0; i<xdimown; i++) {
      z1[i] = x[i] / s;
      z2[i] = y[i] / s;
    }

      //printf("h/s= %10g\n", SQRT((z1[0] - z2[0]) * (z1[0] - z2[0])  + (z1[1] - z2[1]) * (z1[1] - z2[1])));
     //  if (x[0]!=1.0) printf("s=%10g, %10g %10g; %10g %10g \n", s, z1[0], z1[1], z2[0], z2[1]);
    //    assert(z2[1] < 0.3);
	     
  } else if ((notdelete = aniso==NULL && (scale==NULL || scale[0] == 1.0))) {
    z1 = x;
    z2 = y;
  } else {
    dimz = cov->ncol[DANISO];
    assert(xdimown == cov->nrow[DANISO]);
    double *xz1, *xz2;
    TALLOC_GLOBAL_X1(z1, dimz);
    TALLOC_GLOBAL_X2(z2, dimz);
    if (aniso != NULL) {
      xA(x, y, aniso, xdimown, dimz, z1, z2);
      xz1 = z1;
      xz2 = z2;
    } else {
      xz1 = x;
      xz2 = y;
    }
    if (scale != NULL) {
      double s = scale[0];
      if (s > 0.0) {
	double invscale = 1.0 / s;
	for (i=0; i<xdimown; i++) {
	  z1[i] = invscale * xz1[i];
	  z2[i] = invscale * xz2[i];
	}
      } else {
	for (i=0; i<nproj; i++) {
	  z1[i] = (xz1[i] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
	  z2[i] = (xz2[i] == 0.0 && s == 0.0) ? 0.0 : RF_INF;
	}
      }
    }
  }
   

  double var;
  if (cov->Sdollar->simplevar) {
    var = P0(DVAR);
    if (Sign != NULL) var = LOG(var);
  } else {
    assert(equalsCoordinateSystem(OWNISO(0)))
    double w;    
    FCTN(x, info, cov->kappasub[DVAR], &var);
    FCTN(y, infoY, cov->kappasub[DVAR], &w);
    var *= w;
    if (Sign == NULL) var = SQRT(var); else var = 0.5 * LOG(var);    
  }

  if (Scale != NULL) {
    double s12 = s1 * s2 /smeanSq,
      dimhalf = 0.5 * xdimown;
    if (Sign == NULL) 
      if (xdimown == 2) var *= s12; // deal most frequent case faster
      else var *= POW(s12, dimhalf);
    else var += dimhalf * LOG(s12); 
  }

  if (Sign == NULL) {
    NONSTATCOV(z1, z2, info, next, v);
    for (i=0; i<vdimSq; i++) v[i] *= var;
  } else { // Sign != NULL
    LOGNONSTATCOV(z1, z2, info, next, v, Sign);
    for (i=0; i<vdimSq; i++) v[i] += var;
  }

  if (!notdelete) {
    FREE_TALLOC(z1);
    FREE_TALLOC(z2);
  } 
}


void covmatrixS(model *cov, bool ignore_y, double *v) {
  model *next = cov->sub[DOLLAR_SUB];
  assert(next != NULL);
  int 
    dim = Loctsdim(cov),
    vdim = VDIM0;
  assert(dim == PREVTOTALXDIM);

   SPRINTF(cov->base->error_location, "'%.50s'", NICK(cov));

   if ((!PisNULL(DSCALE) && P0(DSCALE) != 1.0) || 
       !PisNULL(DANISO) || !PisNULL(DPROJ) || 
       cov->kappasub[DSCALE] != NULL ||
       cov->kappasub[DAUSER] != NULL ||
       cov->kappasub[DANISO] != NULL ||
       cov->kappasub[DPROJ] != NULL ||
       !DefList[NEXTNR].is_covmatrix(next)
       ) {  
     StandardCovMatrix(cov, ignore_y, v); 
    return;
  }
 
  if (cov->Spgs == NULL && alloc_fctn(cov, dim, vdim * vdim) != NOERROR)
    XERR(ERRORMEMORYALLOCATION);

  int L = OWNLASTSYSTEM; if (L != LASTSYSTEM(PREVSYSOF(next))) BUG;
  for (int s=0; s<=L; s++) {
    if (XDIM(PREVSYSOF(next), s) != NEXTXDIM(s)) {
      BUG; // fuehrt zum richtigen Resultat, sollte aber nicht
      // vorkommen!
      CovarianceMatrix(cov, ignore_y, v); 
      return;
    }
  }

  // assert(locnext != NULL);
  int tot = vdim * LoctotalpointsY(cov, ignore_y),
    totSq = tot * tot;

// save trafos
  Systems_type prev, gatter, own;
  COPYALLSYSTEMS(prev, PREVSYSOF(next), false);
  COPYALLSYSTEMS(gatter, GATTERSYSOF(next), false);
  COPYALLSYSTEMS(own, NEXT, false);

  // next trafos reset by trafos of cov
  COPYALLSYSTEMS(PREVSYSOF(next), PREV, false);
  COPYALLSYSTEMS(GATTERSYSOF(next), GATTER, false);
  COPYALLSYSTEMS(NEXT, OWN, true);
  DefList[NEXTNR].covmatrix(next, ignore_y, v);//hier wird uU next->totalpoints gesetzt

  // restore trafos
  COPYALLSYSTEMS(PREVSYSOF(next), prev, false);
  COPYALLSYSTEMS(GATTERSYSOF(next), gatter, false);
  COPYALLSYSTEMS(NEXT, own, false);
  
  if (cov->Sdollar->simplevar) {
    double var = P0(DVAR);
    if (var == 1.0) return;
    for (int i=0; i<totSq; v[i++] *= var);
  } else {
    BUG;
  }
}

char iscovmatrixS(model *cov) {
  model *sub = cov->sub[DOLLAR_SUB];
  return (int) ((PisNULL(DSCALE) || P0(DSCALE) ==1.0) &&	
		PisNULL( DAUSER) && PisNULL( DANISO) &&
		PisNULL(DPROJ) &&
		cov->Sdollar->simplevar && // to do
		PisNULL(DANISO)) * DefList[SUBNR].is_covmatrix(sub);
}

void DS(double *x, int *info, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert( cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL && cov->kappasub[DSCALE] == NULL);
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj;
  double y[2], varSc,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;
  assert(cov->Sdollar->simplevar);
  assert(isCartesian(OWNISO(0)));

  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varSc = P0(DVAR) * spinvscale;

  if (nproj == 0) {
    y[0] = x[0] * spinvscale; 
    y[1] = (equalsIsotropic(OWNISO(0)) || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    assert(equalsSpaceIsotropic(OWNISO(0)));
    //  assert(P0 INT(DPROJ) == PROJ_SPACE); // #define Nproj
    y[0] = x[0] * spinvscale;
    y[1] = RF_NAN;
  }

  Abl1(y, info, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varSc; 
}

void DDS(double *x, int *info, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL );
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj,
    *proj = PPROJ;
  double y[2], varScSq,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;
  
  assert(isCartesian(OWNISO(0)));
  assert(cov->Sdollar->simplevar);
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScSq = P0(DVAR) * spinvscale * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (equalsIsotropic(OWNISO(0)) || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i]] * x[proj[i]];
    }
    y[0] = SQRT(y[0]) * spinvscale;
  }
  Abl2(y, info, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScSq; 
}


void D3S(double *x, int *info, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL);
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj,
    *proj = PPROJ;
  double y[2], varScS3,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;

  assert(isCartesian(OWNISO(0)));
  assert(cov->Sdollar->simplevar);
  if (aniso != NULL) {
    spinvscale *= aniso[0];
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScS3 = P0(DVAR) * spinvscale * spinvscale * spinvscale;
  
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (OWNISO(0)==ISOTROPIC || cov->ncol[DANISO]==1) ? 0.0
      : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i]] * x[proj[i]];
    }
    y[0] = SQRT(y[0]) * spinvscale;
  }
  Abl3(y, info, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScS3; 
}

void D4S(double *x, int *info, model *cov, double *v){
  model *next = cov->sub[DOLLAR_SUB];
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL);
  int i,
    vdim = VDIM0,
    vdimSq = vdim * vdim,
    nproj = Nproj,
    *proj = PPROJ;
  double y[2], varScS4,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    spinvscale = 1.0;

  assert(isCartesian(OWNISO(0)));
  if (aniso != NULL) {
    spinvscale *= aniso[0]; 
    // was passiert, wenn es aniso nicht vom TypeIso ist ??
  }
  if (scale != NULL) spinvscale /= scale[0];
  varScS4 = spinvscale * spinvscale;
  varScS4 *= varScS4 * P0(DVAR);
  if (nproj == 0) {
    y[0] = x[0] * spinvscale;
    y[1] = (equalsIsotropic(OWNISO(0)) || cov->ncol[DANISO]==1) ? 0.0
     : x[1] * aniso[3]; // temporal; temporal scale
  } else {
    BUG;
    for (i=0; i<nproj; i++) {
      y[0] += x[proj[i]] * x[proj[i]];
    }
    y[0] = SQRT(y[0]) * spinvscale;
  }
  Abl4(y, info, next, v);
  for (i=0; i<vdimSq; i++) v[i] *= varScS4; 
}


void nablahessS(double *x, int *info, model *cov, double *v, bool nabla){
  model *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER] :
    cov->kappasub[DANISO];
  if (Aniso != NULL) BUG;
  int i, endfor,
    dim = cov->nrow[DANISO],// == ncol == x d i m ?
    xdimown = OWNTOTALXDIM,
    nproj = Nproj;
  double *xy, *vw,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    var = P0(DVAR);
  if (nproj != 0) BUG;
  if (dim != xdimown) BUG;
	   
  if (!cov->Sdollar->simplevar) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");


  TALLOC_DOUBLE(z0);
  TALLOC_DOUBLE(w0);
  if (aniso != NULL) {  
    TALLOC_GLOBAL_X1(z0, xdimown);
    TALLOC_GLOBAL_X2(w0, xdimown);
    xA(x, aniso, xdimown, xdimown, z0);
    xy = z0;
    vw = w0;
  } else {
    xy = x;
    vw = v;
  }

  TALLOC_DOUBLE(y0);
  if (scale != NULL) {
    TALLOC_GLOBAL_X3(y0, xdimown);
    assert(scale[0] > 0.0);
    double spinvscale = 1.0 / scale[0];
    var *= spinvscale;
    if (!nabla) var *= spinvscale; // gives spinvscale^2 in case of hess
    for (i=0; i<xdimown; i++) y0[i] = xy[i] * spinvscale;
    xy = y0;
  }

  endfor = xdimown;
  if (nabla) {
    NABLA(xy, info, next, vw);
  } else {
    HESSE(xy, info, next, vw);
    endfor *= xdimown;
  }
     
  if (aniso != NULL) {  
    if (nabla) Ax(aniso, vw, xdimown, xdimown, v); 
    else {
      XCXt(aniso, vw, v, xdimown, xdimown); //todo:?reprogramm XCXt with allocation here ?
    }
    FREE_TALLOC(z0);
    FREE_TALLOC(w0);
  }

  FREE_TALLOC(y0);
  
  for (i=0; i<endfor; i++) v[i] *= var; 
}

void nablaS(double *x, int *info, model *cov, double *v){
  nablahessS(x, info, cov, v, true);
}
void hessS(double *x, int *info, model *cov, double *v){
  nablahessS(x, info, cov, v, false);
}


 

void inverseS(double *x, model *cov, double *v) {
  model *next = cov->sub[DOLLAR_SUB];
  int i,
    idx[3] = {DANISO, DAUSER, DPROJ};
  double 
    scale;
  
  if (cov->kappasub[DVAR] != NULL) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");

  for (i=0; i<3; i++) {
    if (cov->kappasub[idx[i]] != NULL) {
      ERR1( "inverse can only be calculated if '%.20s' is not an arbitrary function",
	    KNAME(idx[i])); 
    }
  }
  if (cov->kappasub[DSCALE] != NULL) {
    double left;
    model *sub = cov->kappasub[DSCALE];
    INVERSENONSTAT(ZERO(sub), sub, &left, &scale);
     if (left < 0.0) ERR("scale negative.");
  } else scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
 
  int
    dim = PREVTOTALXDIM,
    nproj = Nproj;
  //    *proj = (int *)P(DPROJ];
  double y, 
    s = 1.0,
    *aniso=P(DANISO),
    var = P0(DVAR);

  if (dim != 1) BUG;

  if (aniso != NULL) {
    if (isMiso(Type(aniso, cov->nrow[DANISO], cov->ncol[DANISO]))) 
      s /= aniso[0];
    else NotProgrammedYet(""); // to do
  }
  s *= scale;  
  if (nproj == 0) {
    y= *x / var; // inversion, so variance becomes scale
  } else {
    BUG;  //ERR("nproj is not allowed in invS");
  }
  
  if (DefList[NEXTNR].inverse == inverseErr) BUG;
  INVERSE(&y, next, v);
 
  for (i=0; i<dim; i++) v[i] *= s; //!

}


void inversenonstatS(double *x, model *cov, double *left, double*right,
		     bool log){
  model
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE];
   int 
    dim = PREVTOTALXDIM,
    nproj = Nproj;
  //    *proj = (int *)P(DPROJ];
  double y, 
    s = 1.0,
    *scale =P(DSCALE),
    *aniso=P(DANISO),
    var = P0(DVAR);
 
  if (cov->kappasub[DVAR] != NULL) 
    NotProgrammedYet("nabla not programmed for arbitrary 'var'");


  if (nproj == 0) {
    y= *x / var; // inversion, so variance becomes scale
  } else {
    BUG;  //ERR("nproj is not allowed in invS");
  }

  
  if (DefList[NEXTNR].inverse_nonstat == inversenonstatErr) BUG;
  if (log) {
    INVERSENONSTATLOG(&y, next, left, right);
    //   printf("logS %f  %f %f %s\n", y, *left, *right, NAME(next));    
  } else {
    INVERSENONSTAT(&y, next, left, right);
    //    printf("S %f  %e %e x=%f, var=%f  %s\n", y, *left, *right, *x, var, NAME(next));
    //if (y > 1.0) {PMI(cov); crash();}
  }
  
  if (aniso != NULL) {
    if (isMiso(Type(aniso, cov->nrow[DANISO], cov->ncol[DANISO]))) s/=aniso[0];
    else {
      getStorage(S ,       dollar);
      int   
	ncol = cov->ncol[DANISO],
	nrow = cov->nrow[DANISO],
	ncnr = ncol * nrow,
	size = ncol * sizeof(double);
      bool redo;
      Long bytes = (Long) ncnr * sizeof(double);
      if (ncol != nrow) BUG;
      redo = S->save_aniso == NULL;
      ALLC_NEW(Sdollar, save_aniso, ncnr, save_aniso);
      ALLC_NEW(Sdollar, inv_aniso, ncnr, inv_aniso);
      TALLOC_X1(LR, ncol);
      S->done = false;
      int pid;
      if (!redo) {
	for (int i=0; i<ncnr;i++)
	  if ((redo = save_aniso[i] != P(DANISO)[i])) break;
	Ext_pid(&pid);
	if (S->pid == 0) {
	  S->pid = pid;	  
	  if (S->pid == pid) { // two processes writing on the same place
	    S->busy = true; // one process doing the job
	  } else if (cov->base->global_utils.basic.cores == 1 || !S->busy)
	    ERR("Multiprocessor conflict. Plse inform maintainer & try again");
	} else if (cov->base->global_utils.basic.cores == 1 || !S->busy) 
	  ERR("Multiprocessor conflict. Plse inform maintainer & try again");
      }
      if (redo && S->pid == pid) {
	MEMCOPY(save_aniso, P(DANISO), bytes);
	MEMCOPY(inv_aniso, P(DANISO), bytes);
	if (Ext_invertMatrix(inv_aniso, nrow) != NOERROR) XERR(ERRORANISO_INV);
	S->busy = false;
	S->done = true;
      }
      while(S->busy && !S->done) {
	int sleep = (double) nrow * nrow * nrow / 1e8;
	sleep = sleep == 0 ? 1 : sleep;
	if (PL >= PL_BRANCHING) { PRINTF("|");}
	Ext_sleepMicro(&sleep);
      }
      
      MEMCOPY(LR, right, size);
      xA(LR, inv_aniso, nrow, ncol, right);
      
      MEMCOPY(LR, left, size);
      xA(LR, inv_aniso, nrow, ncol, left);      
      END_TALLOC_X1;
    }
  }

  if (Aniso != NULL) {
    if (aniso != NULL) BUG;
    if (DefList[MODELNR(Aniso)].inverse == inverseErr) XERR(ERRORANISO_INV);
    int 
      nrow = Aniso->vdim[0],
      ncol = Aniso->vdim[1],
      size = nrow * sizeof(double);
    if (dim != ncol || ncol != 1)
      ERR("anisotropy function not of appropriate form");
    TALLOC_X1(LR, nrow);
    
    MEMCOPY(LR, right, size);
    INVERSE(LR, Aniso, right);
          
    MEMCOPY(LR, left, size);
    INVERSE(LR, Aniso, left);
    END_TALLOC_X1;
  }

  if (Scale != NULL && !isnowRandom(Scale)) {
    double dummy;
    Zero(Scale, &dummy);
    s *= dummy;
  } else {
    if (scale != NULL) s *= scale[0];

    // printf("name %s\n", DefList[NEXTNR].name);
    int nextdim = NEXTTOTALXDIM;
    if (dim > nextdim && isAnyIsotropic(NEXTISO(0)) &&
	aniso == NULL && Aniso == NULL) {
      assert(nextdim == 1);
      double r = right[0],
	l = left[0];
      //   printf("%e %e\n", l, r);
      if (-l != r) BUG;
      for (int i=1; i<dim; i++) {
	left[i] = l;
	right[i] = r;
      }
    }
  }
  if (s != 1.0) {
    //    PMI0(cov);
    //printf("dim = %d\n", dim);
    //BUG;
    for (int i=0; i<dim; i++) {
      left[i] *= s; //!
      right[i] *= s;
    }
  }

  //  PMI0(cov);
  //  printf("SS %f  %e %e %f dim=%d %s\n", y, *left, right[0], s, dim, NAME(next));
 
}

void inversenonstatS(double *x, model *cov, double *left, double*right) {
  inversenonstatS(x, cov, left, right, false);
}

void inverse_log_nonstatS(double *x, model *cov, double *left, double*right){
 inversenonstatS(x, cov, left, right, true);
}

void coinitS(model *cov, localfactstype *li) {
  assert(cov->Sdollar->simplevar);
  model *next = cov->sub[DOLLAR_SUB];
  if ( DefList[NEXTNR].coinit == NULL)
    ERR("# cannot find coinit -- please inform author");
  DefList[NEXTNR].coinit(next, li);
}
void ieinitS(model *cov, localfactstype *li) {
   assert(cov->Sdollar->simplevar);
 model *next = cov->sub[DOLLAR_SUB];
  
  if ( DefList[NEXTNR].ieinit == NULL)
    ERR("# cannot find ieinit -- please inform author");
  DefList[NEXTNR].ieinit(next, li);
}

void tbm2S(double *x, int *info, model *cov, double *v){
  assert(cov->Sdollar->simplevar);
 model *next = cov->sub[DOLLAR_SUB];
 // defn *C = DefList + NEXTNR; // kein gatternr, da isotrop
  double y[2],  *xy,
    *scale =P(DSCALE),
    *aniso = P(DANISO);
  assert(cov->kappasub[DAUSER] == NULL && cov->kappasub[DANISO] == NULL);
 
  assert(Nproj == 0);
  if (aniso!=NULL) {
    if (cov->ncol[DANISO]==2) {  // turning layers
      y[0] = x[0] * aniso[0]; // spatial 
      y[1] = x[1] * aniso[3]; // temporal
      assert(aniso[1] == 0.0 && aniso[2] == 0.0);
      if (y[0] == 0.0) COV(y, info, next, v) else TBM2CALL(y, info, next, v);
    } else {
      assert(cov->ncol[DANISO]==1);
      if (cov->nrow[DANISO] == 1) { // turning bands
	y[0] = x[0] * aniso[0]; // purely spatial 
 	TBM2CALL(y, info, next, v);
      } else { // turning layers, dimension reduction
	if (P0(DANISO) == 0.0) {
	  y[0] = x[1] * aniso[1]; // temporal 
	  COV(y, info, next, v); 
	} else {
	  y[0] = x[0] * aniso[0]; // temporal 
	  TBM2CALL(y, info, next, v);
	}
      } 
    }
    xy = y;
  } else xy = x;

   if (scale != NULL) {
    double s = scale[0];
    if (s > 0) { 
      double invscale = 1.0 / s;
      if (OWNTOTALXDIM == 2){
	y[0] = xy[0] * invscale; // spatial 
	y[1] = xy[1] * invscale; // temporal
	if (y[0] == 0.0) COV(y, info, next, v) else TBM2CALL(y, info, next, v);
      } else {
	y[0] = xy[0] * invscale; // purely spatial 
	TBM2CALL(y, info, next, v);
     }
   } else {
      y[0] = (s < 0.0 || xy[0] != 0.0) ? RF_INF : 0.0;
      if (OWNTOTALXDIM == 2)
 	y[1] = (s < 0.0 || xy[1] != 0.0) ? RF_INF : 0.0;
      TBM2CALL(y, info, next, v);
    }
   } else {
     TBM2CALL(xy, info, next, v);
   }
  *v *= P0(DVAR);
}


// TODO : Aniso=Matrix: direkte implementierung in S,
// sodass nicht ueber initS gegangen werden muss, sondern
//  e  < -  e^\top Aniso


int TaylorS(model *cov) {
  model 
    *next = cov->sub[DOLLAR_SUB],
    *sub = cov->key == NULL ? next : cov->key;
  int i;

  if (PisNULL(DPROJ) && PisNULL(DANISO)) {
    double scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
    cov->taylorN = sub->taylorN;  
    for (i=0; i<cov->taylorN; i++) {
      cov->taylor[i][TaylorPow] = sub->taylor[i][TaylorPow];
      cov->taylor[i][TaylorConst] = sub->taylor[i][TaylorConst] *
	P0(DVAR) * POW(scale, -sub->taylor[i][TaylorPow]);   
    }

    cov->tailN = sub->tailN;  
    for (i=0; i<cov->tailN; i++) {
      cov->tail[i][TaylorPow] = sub->tail[i][TaylorPow];
      cov->tail[i][TaylorExpPow] = sub->tail[i][TaylorExpPow];
      cov->tail[i][TaylorConst] = sub->tail[i][TaylorConst] *
	P0(DVAR) * POW(scale, -sub->tail[i][TaylorPow]);   
      cov->tail[i][TaylorExpConst] = sub->tail[i][TaylorExpConst] *
	POW(scale, -sub->tail[i][TaylorExpPow]);
    }
  } else {
    cov->taylorN = cov->tailN = 0;
  }
  RETURN_NOERROR;
}

int checkS(model *cov) {
  globalparam *global = &(cov->base->global);
  //  printf("checkS!!\n");

  // hier kommt unerwartet  ein scale == nan rein ?!!
  model 
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE],
    *Var = cov->kappasub[DVAR],
    *sub = cov->key == NULL ? next : cov->key;
   isotropy_type owniso = OWNISO(0);
   int  err,
     //xdimgatter = GATTERTOTALXDIM,
     // KOMMENTAR NICHT LOESCHEN
     // *proj = PPROJ, // auf keinen Fall setzen, da Pointer unten neu
   //                        gesetzt wird!!!!
     last = OWNLASTSYSTEM; 
  matrix_type mtype = TypeMany;
  bool isProcess = COVNR == DOLLAR_PROC;

  // if (Aniso == NULL && Scale == NULL && Var == NULL &&
  //  PisNULL(DANISO) &&  PisNULL(DAUSER) &&  PisNULL(DSCALE) &&  PisNULL(DVAR))
  //  SERR("no parameter is set.");
 
  assert(isAnyDollar(cov));
  if (!isDollarProc(cov)) set_nr(OWN, FIRSTDOLLAR); //wegen nr++ unten! NICHT SET_NR!
  
  // cov->q[1] not NULL then A has been given
  
  if ((err = checkkappas(cov, false)) != NOERROR) { // muss weit vorne stehen!
    RETURN_ERR(err);
  }
  kdefault(cov, DVAR, 1.0);
 //  printf("checkS2!!\n");

  if ( (!PisNULL(DSCALE) || Scale != NULL) && OWNLASTSYSTEM != 0)
    SERR1("'%.50s' may be used only if a single coordinate system is envolved", 
	  KNAME(DSCALE));

  bool angle = isAngle(Aniso);
  int dim = OWNLOGDIM(0);

  //ownxdim bei allen processen kritisch suchen: ist es hier ok#>>>

  assert(!isProcess || dim == OWNXDIM(0));
  
  if (angle && PisNULL(DANISO) && PisNULL(DAUSER) &&
      OWNLASTSYSTEM == 0 && OWNXDIM(0) == dim) {
    ASSERT_CARTESIAN;
    if (!isCartesian(owniso)) RETURN_ERR(ERRORANISO);
    if ((err = CHECK(Aniso, dim, dim, ShapeType, XONLY, CARTESIAN_COORD,
		     SUBMODEL_DEP, cov->frame)) == NOERROR) {
      if (isverysimple(Aniso)) {     
 	 PALLOC(DAUSER, dim, dim);
 	 AngleMatrix(Aniso, P(DAUSER));
	 COV_DELETE(cov->kappasub + DAUSER, cov);
	 Aniso = cov->kappasub[DAUSER] = NULL;
	 angle = false;
      }
    }
  }

  if (cov->kappasub[DANISO] != NULL) // kann auch weggelassen werden. sollte aufgefangen werden
    SERR1("'%.50s' may not be a function.", KNAME(DANISO));
  if (!PisNULL(DAUSER)) {
    if (!isCartesian(owniso))  RETURN_ERR(ERRORANISO);
    //    if (GLO BAL.messages.warn_Aniso) {
    //      PRINTF("NOTE! Starting with RandomFields 3.0, the use of '%s' is different from\nthe former '%s' insofar that '%s' is multiplied from the right by 'x' (i.e. Ax),\nwhereas '%s' had been multiplied from the left by 'x' (i.e. xA).\n", KNAME(DAUSER), KNAME(DANISO), KNAME(DANISO), KNAME(DAUSER));
    //    }  GLO BAL.messages.warn_Aniso = false;
    
    // here for the first time
    if (!PisNULL(DANISO)) RETURN_ERR(ERRORANISO_T); 
    int 
      lnrow = cov->nrow[DAUSER],
      lncol = cov->ncol[DAUSER];
    Long
      total = lncol * lnrow;
	
    double
      *pA = P(DAUSER); 
    PALLOC(DANISO, lncol, lnrow); // !! ACHTUNG col, row gekreuzt
    for (int i=0, k=0; i<lnrow; i++) {
      for (Long j=i; j<total; j+=lnrow) {
	if (ISNAN(pA[j])) ERR("(Parts of) 'Aniso' may not be estimated. Instead 'anisoT', which equals 't(Aniso)', must be used. See ?RMS for details.");
	P(DANISO)[k++] = pA[j];
      }
    }
    PFREE(DAUSER);
  }

  bool simplevar = Var == NULL || isnowRandom(Var); // checkkappa s.o.
  if (!simplevar) {
    ptwise_type ptt = cov->ptwise_definite;
    assert(equalsCoordinateSystem(owniso));
    if ((err = CHECK(Var,
		     OWNLOGDIM(0), OWNXDIM(0), 
		     ShapeType, // only!! -- for pos def use RMprod
		     XONLY, owniso,
		     SCALAR, // s.a. ptwise_posdef unten!
		     TrendType)) != NOERROR) RETURN_ERR(err);
    if (Var->ptwise_definite != pt_posdef) {
      if (Var->ptwise_definite == pt_unknown) {
	if (global->messages.warn_negvar && cov->q==NULL) { // just a flag
	  QALLOC(1);
	  PRINTF("Positivity of the variance in '%.50s' cannot be detected.\n",
		  NICK(next));
	}
      } else {
	SERR2("positivity of '%.50s' required. Got '%.50s'", KNAME(DVAR), 
	      POSITIVITY_NAMES[Var->ptwise_definite]);
      }
    }
    cov->ptwise_definite = ptt;
  }

  ONCE_NEW_STORAGE(dollar);
  ONCE_EXTRA_STORAGE;

  int xdimNeu = OWNXDIM(0);
 //  printf("checkS ff !!\n");
  if (Aniso != NULL) { // CASE 1: Aniso is given by arbitrary fctn
    if (!isDollarProc(cov) && !angle && !isKernel(OWN)) RETURN_ERR(ERRORFAILED);
    if (!PisNULL(DANISO) || !PisNULL(DPROJ) || !PisNULL(DSCALE) ||
	(Scale!=NULL && !isnowRandom(Scale)) )
      SERR2("if '%.50s' is an arbitrary function, only '%.50s' may be given as additional parameter.", KNAME(DAUSER), KNAME(DVAR));

    //assert(equalsCartCoord(owniso) || (angle && equalsSymmetric(owniso)));
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = wahr;
    if ((err = CHECK(Aniso, OWNLOGDIM(0), OWNXDIM(0), ShapeType, XONLY,
		     CARTESIAN_COORD, SUBMODEL_DEP, cov->frame)) != NOERROR) {
       RETURN_ERR(err);
    }

    if (Aniso->vdim[1] != 1)
      SERR3("'%.50s' returns a matrix of dimension %d x %d, but a vector is need.",
	    KNAME(DANISO), Aniso->vdim[0], Aniso->vdim[1]);

    if (cov->key==NULL) {
      ASSERT_QUASIONESYSTEM;
      if ((err = CHECK(sub, Aniso->vdim[0], Aniso->vdim[0], 
		      OWNTYPE(0), OWNDOM(0), owniso, 
		      SUBMODEL_DEP, cov->frame)) != NOERROR) {
	RETURN_ERR(err);
      }
    }
    cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = 
      cov->pref[Hyperplane] = cov->pref[SpectralTBM] = cov->pref[TBM] = 
      sub->pref[CircEmbed] = sub->pref[CircEmbedCutoff] = 
      sub->pref[CircEmbedIntrinsic] = sub->pref[Sequential] = 
      sub->pref[Specific] = PREF_NONE;

  } else if (Scale != NULL && !isnowRandom(Scale)) {
     // CASE 2: Scale is given by arbitrary function
    if (!isDollarProc(cov) && !isKernel(OWN)) RETURN_ERR(ERRORFAILED);
    if (!PisNULL(DANISO) || !PisNULL(DPROJ) || !PisNULL(DSCALE))
       SERR2("if '%.50s' is an arbitrary function, only '%.50s' may be given as additional parameter.", KNAME(DSCALE), KNAME(DVAR));
 
    assert(equalsCartCoord(owniso));
    
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = wahr;
    assert(equalsCartCoord(owniso));

    if ((err = CHECK(Scale, OWNLOGDIM(0), OWNXDIM(0), ShapeType, XONLY,
		     owniso, SCALAR, TrendType)) != NOERROR) {
      RETURN_ERR(err);
    }
    if (Scale->vdim[1] != 1 || Scale->vdim[0] != 1)
      SERR3("'%.50s' must be scalar, not %d x %d",
	    KNAME(DSCALE), Aniso->vdim[0], Aniso->vdim[1]);
    if (cov->key==NULL) {
      ASSERT_QUASIONESYSTEM;
      if ((err = CHECK(sub, OWNLOGDIM(0), OWNXDIM(0),
		       OWNTYPE(0), OWNDOM(0), 
		       owniso, SUBMODEL_DEP, cov->frame)) != NOERROR) {
	RETURN_ERR(err);
      }
      if (!isNormalMixture(next)) 
	SERR("scale function only allowed for normal mixtures.");
    }
    
 
    cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = 
      cov->pref[Hyperplane] = cov->pref[SpectralTBM] = cov->pref[TBM] = 
      sub->pref[CircEmbed] = sub->pref[CircEmbedCutoff] = 
      sub->pref[CircEmbedIntrinsic] = sub->pref[Sequential] = 
      sub->pref[Specific] = PREF_NONE;
    
  } else if (!PisNULL(DANISO)) { // CASE 3, anisotropy matrix is given
    int 
      nrow = cov->nrow[DANISO],
      ncol = cov->ncol[DANISO];

   
    if (nrow==0 || ncol==0) SERR("dimension of the matrix is 0");
    if (!PisNULL(DPROJ)) RETURN_ERR(ERRORANISO_T);

    int xdimown = OWNTOTALXDIM;

    if (xdimown != nrow) {
      if (PL >= PL_ERRORS) {LPRINT("xdim=%d != nrow=%d\n", xdimown, nrow);}
      SERR("#rows of anisotropy matrix does not match dim. of coordinates");
    }
    if (xdimown != OWNLOGDIM(0) && nrow != ncol)
      SERR("non-quadratic anisotropy matrices only allowed if dimension of coordinates equals spatio-temporal dimension");

    mtype = Type(P(DANISO), nrow, ncol);
    
    cov->full_derivs = cov->rese_derivs = 0;
    cov->loggiven = wahr;

    if (last == 0 || isSpaceIsotropic(OWN)) {
      if (!equalsSpaceIsotropic(OWN)) {
	switch (owniso) {
	case ISOTROPIC :   
	  if (OWNXDIM(0) != 1) RETURN_ERR(ERRORANISO);
	  cov->full_derivs = cov->rese_derivs = 2;
	  break;
	case VECTORISOTROPIC :
	  if (!isMiso(mtype)) RETURN_ERR(ERRORANISO); 
	break;
	case SYMMETRIC: case CARTESIAN_COORD :
	  break;
	case PREVMODEL_I : BUG;      
	case GNOMONIC_PROJ :  case ORTHOGRAPHIC_PROJ:
	  if (!isnowProcess(cov)) RETURN_ERR(ERRORANISO);
	  break;
	default :
	  if (isCartesian(owniso)) { BUG; }
	  else RETURN_ERR(ERRORANISO);
	}
	
	if (!isMiso(mtype)) cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;
	ASSERT_ONESYSTEM;
	err = CHECK(sub, ncol, ncol, OWNTYPE(0), OWNDOM(0), 
		    ncol==1 && !isnowProcess(cov)
		    ? IsotropicOf(owniso) : owniso, 
		    SUBMODEL_DEP, cov->frame);
      } else {	
	// spaceisotropic
	/*
    
	  case DOUBLEISOTROPIC :  
	  cov->full_derivs =  cov->rese_derivs = isMdiag(mtype) ? 2 : 0;
	  if (nrow != 2 || !isMdiag(mtype)) {
	  SERR("spaceisotropy needs a 2x2 diagonal matrix");
	  }
	  break;    	  
	*/
	
	cov->full_derivs = cov->rese_derivs = isMdiag(mtype) ? 2 : 0;
	if (nrow != 2 || !isMdiag(mtype)) 
	  SERR("spaceisotropy needs a 2x2 diagonal matrix");
	if (nrow != ncol) BUG;
	COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
	err = CHECK_GEN(sub, SUBMODEL_DEP, SUBMODEL_DEP, cov->frame, false);
      }

      if (err != NOERROR) RETURN_ERR(err);

    } else { // more than 1 coordinate system
      if (!isMdiag(mtype) ) 
	SERR1("If several coordinate systems are envolved, currently '%.50s' can only be a diagonal matrix", KNAME(DANISO));
      
      COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
      if ((err = CHECK_GEN(sub, SUBMODEL_DEP, SUBMODEL_DEP, cov->frame, false)) 
	  != NOERROR) {
	RETURN_ERR(err);
      }
      cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;

  /*

    if (!isMdiag(mtype)) 
      cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] = cov->pref[Hyperplane] = PREF_NONE;
    if (owniso != DOUBLEISOTROPIC && !isMiso(mtype))
      cov->pref[SpectralTBM] = cov->pref[TBM] = PREF_NONE;

      if ((err = CHECK(suib, ncol, ncol, cov->typus, OWNDOM(0), 
		 ncol==1 && !isProcess(cov->typus) ? ISOTROPIC : owniso, 
		 SUBMODEL_DEP, cov->frame))
	  != NOERROR) {
 	RETURN_ERR(err);
      }
    */
      
    }
 
    if (!isMdiag(mtype)) 
      cov->pref[Nugget] = cov->pref[RandomCoin] = cov->pref[Average] =
	cov->pref[Hyperplane] = PREF_NONE;
  
  } else if (!PisNULL(DPROJ)) { // geaendert werden: xdim, logdim, iso
    location_type *loc = LocPrev(cov);
    COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);
    int xdim = LocLocxdimOZ(loc) + LocLocTime(loc),
      *p = PINT(DPROJ);  // OK

    FREE(PPROJ);
    if (p[0] > 0) {
      Nproj = cov->nrow[DPROJ]; // OK
      int bytes = sizeof(int) * Nproj; 
      PPROJ = (int*) MALLOC(bytes);
      for (int i=0; i<Nproj; i++) PPROJ[i] = p[i] - 1;
    } else {
      int nproj = cov->nrow[DPROJ]; // OK
      if (nproj != 1)
	SERR1("values 'space' and 'time' in argument '%.50s' do not allow additional values", KNAME(DPROJ));
      if (!(LocLocTime(loc) && xdim >= 2 &&
	  (isCartesian(OWNISO(last))
#if MAXSYSTEMS == 1
	   || (isAnySpherical(PREVISO(0)) && PREVXDIM(0) > 2)
#endif							      
	   ))) {
	//	PMI0(cov->calling);
	// PMI0(cov);
	// printf("%d %d %d %d\n",  LocLocTime(loc), xdim,   isCartesian(OWNISO(last)),	       isAnySpherical(PREVISO(0)) && PREVXDIM(0) > 2);
	SERR1("unallowed use of '%.50s' or model too complicated.",
	      KNAME(DPROJ));
      }

      bool Space = p[0] == PROJ_SPACE;
      assert(Space || p[0]  == PROJ_TIME);
      Nproj = equalsSpaceIsotropic(OWN) || !Space ? 1 : xdim - 1;
      PPROJ = (int*) MALLOC(sizeof(int) * Nproj);
      assert(!equalsSpaceIsotropic(OWN) || !isProcess);
      if (equalsSpaceIsotropic(OWN)) PPROJ[0] = !Space;
      else if (Space) for (int i=0; i<Nproj; i++) PPROJ[i] = i;
      else PPROJ[0] = xdim - 1; // R coding 

    //      printf("Space %d %ld %d %d\n", Space, p, PROJ_TIME, xdim);
    }

    xdimNeu = Nproj;
    int max = xdim;
    if (xdimNeu == 0) SERR("Projection to null space not allowed.");

    if (xdimNeu < xdim) {	
      cov->pref[CircEmbed] = cov->pref[CircEmbedCutoff] =
	cov->pref[CircEmbedIntrinsic] =  cov->pref[Direct] =
	cov->pref[Sequential] = 1;
      cov->pref[TBM] = cov->pref[SpectralTBM] = cov->pref[Average]
	= cov->pref[Nugget]= cov->pref[RandomCoin] = cov->pref[Hyperplane]
	=  PREF_NONE;
      cov->pref[Specific] = PREF_BEST;
    }
     
    /// hier Aenderung ab version 3.2.2 so dass spaceisotropic
    // nicht mehr ueber 2 coordinatensysteme dargestellt werden kann
    int nproj = Nproj;
    for (int i=0; i<nproj; i++) {
      int idx = PPROJ[i];
      if (idx < 0)
	SERR1("only positive values allowed for '%.50s'", KNAME(DPROJ))
	else if (idx >= max) SERR2("value of %.50s[%d] too large", KNAME(DPROJ),i);
	
      for (int j=i+1; j<nproj; j++) {
	int jdx = PPROJ[j];
	if (jdx == idx) SERR1("values of '%.50s' must be distinct.", KNAME(DPROJ));
      }
    }


    set_xdim(PREVSYSOF(sub), 0, xdimNeu);
    set_logdim(PREVSYSOF(sub), 0, xdimNeu);
    //    PMI(cov);    printf("PP = %d %s\n", PPROJ[0], NAME(cov));
    if (equalsSpaceIsotropic(OWN)) {
      // printf("%d %d\n", xdimNeu, nproj);
      if (xdimNeu > 2) SERR("maximum length of projection vector is 2");
      if (xdimNeu == 2) {
	if (PPROJ[0] >= PPROJ[1])
	  SERR1("in case of '%.50s' projection directions must be ordered",
		ISO_NAMES[DOUBLEISOTROPIC]);
      } else { 
	assert(xdimNeu == 1);
#if MAXSSYSTEMS == 1
	RETURN_ERR(ERRORWRONGISO);
#else	    
	set_iso(PREVSYSOF(sub), 0, ISOTROPIC);
	if (PPROJ[0] == 0) { // space
	  set_logdim(PREVSYSOF(sub), 0, OWNLOGDIM(0) - 1);
	} else { // time
	  assert(PPROJ[0] == 1);
	  set_logdim(PREVSYSOF(sub), 0, 1);
	}
#endif	    
      }
    } else {// ! spaceisotropic
      switch (OWNISO(0)) {
      case ISOTROPIC : case EARTH_ISOTROPIC : case SPHERICAL_ISOTROPIC : 
	if (OWNXDIM(0) != 1) RETURN_ERR(ERRORANISO);
	break;
      case VECTORISOTROPIC :
	SERR("projection of vectorisotropic fields not programmed yet"); // to do: vdim muss auch reduziert werden ... --- das wird
	// grausam !
	break;
      case SYMMETRIC: case CARTESIAN_COORD:
	break;
      case GNOMONIC_PROJ :  case ORTHOGRAPHIC_PROJ :
	set_iso(PREVSYSOF(sub), 0, CARTESIAN_COORD);
	break;
      case PREVMODEL_I : BUG;      
	break;
      case SPHERICAL_SYMMETRIC : case EARTH_SYMMETRIC :
	if (nproj != 2 || PPROJ[0] != 0 || PPROJ[1]  != 1) {
	  for (int ii=0; ii<nproj; ii++)
	    if (PPROJ[ii] <= 1) APMI0(cov); // RETURN_ERR(ERRORANISO);
	  /// ehemals  owniso = SYMMETRIC; -- ueberall owniso
	  set_iso(PREVSYSOF(sub), 0, SYMMETRIC);
	}
	break;
      case SPHERICAL_COORD : case EARTH_COORD :
	if (nproj != 2 || PPROJ[0] != 0 || PPROJ[1]  !=  1) {
	  for (int ii=0; ii<nproj; ii++) // ??
	    if (PPROJ[ii] <= 1) RETURN_ERR(ERRORANISO);
	  set_iso(PREVSYSOF(sub), 0, CARTESIAN_COORD);
	}
	break;
      default : 
	if (isCartesian(OWNISO(0))) {BUG;}
	else return  ERRORANISO;  // todo
      }
    }
    
      
    //    if ((err = CHECK(sub, nproj, nproj, OWNTYPE(0), OWNDOM(0), owniso,
    //		     SUBMODEL_DEP, cov->frame)) != NOERROR)  RETURN_ERR(err); 
    if ((err = CHECK_GEN(sub, 
			 VDIM0, // SUBMODEL_DEP; geaendert 20.7.14
			 VDIM1, cov->frame, true)) != NOERROR) {
      RETURN_ERR(err);
    }
   
    
  } else { // nproj == 0
 //  printf("checkS  u!!\n");
    COPYALLSYSTEMS(PREVSYSOF(sub), OWN, false);

    // verhindern, dass die coordinaten-transformation anlaeuft,
    // da aus z.B. Erd-coordinaten durch Projektion (auf die Zeitachse)
    // kartesische werden

    
    if ((err = CHECK_GEN(sub, 
			 VDIM0, // SUBMODEL_DEP; geaendert 20.7.14
			 VDIM1, cov->frame, true)) != NOERROR) {
      RETURN_ERR(err);
    }
    
      
    /*   if (cov->key==NULL) {
	 if ((err = CHECK_NO_TRAFO(next, tsdim, xdim Neu, cov->typus, OWNDOM(0),
	 owniso, 
	 VDIM0, // SUBMODEL_DEP; geaendert 20.7.14
	 cov->frame)) != NOERROR) {
	 RETURN_ERR(err);
	 }

	 if (next->domown == OWNDOM(0) &&
	 next->owniso == owniso) // darf nicht nach C-HECK als allgemeine Regel ruebergezogen werden, u.a. nicht wegen stat/nicht-stat wechsel !!
	 // hier geht es, da domown und owniso nur durchgegeben werden und die Werte      // bereits ein Schritt weiter oben richtig/minimal gesetzt werden.
	 } else {
	 if ((err = CHECK_NO_TRAFO(cov->key, tsdim, xdimN eu, cov->typus,
	 OWNDOM(0), owniso,
	 SUBMODEL_DEP, cov->frame)) != NOERROR) RETURN_ERR(err);
	 }

	 */

  } // end no aniso ((end of CASE 4))
//  printf("checkS 664!!\n");
 
  if (( err = checkkappas(cov, false)) != NOERROR) {
    RETURN_ERR(err);
  }
  
 //  printf("checkS 4!!\n");
  setbackward(cov, sub);
  for (int s=0; s<=last; s++) set_maxdim(OWN, s, OWNLOGDIM(s));

  if ((Aniso != NULL || (Scale != NULL && !isnowRandom(Scale)) || 
       !PisNULL(DANISO) || !PisNULL(DPROJ))) { 
    for (int s=0; s<=last; s++) {
      int max = MAXDIM(SUB, s);
      if (MAXDIM(OWN, s) < max) set_maxdim(OWN, s, max);
    }
  }

  if ((!isAnyIsotropic(owniso) && !isDollarProc(cov)) || last > 0){
    // multivariate kann auch xdimNeu == 1 problematisch sein
    set_nr(OWN, COVNR + 1); // JA NICHT SET_NR!!
  }
  // if (!isAnyIsotropic(owniso) && !isDollarProc(cov)) { // multivariate kann auch xdimNeu == 1 problematisch sein
  //   COVNR++;
  // }
  
  if (xdimNeu > 1) {
    cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 0;
  }

  // 30.10.11 kommentiert:
  //  cov->pref[CircEmbedCutoff] = cov->pref[CircEmbedIntrinsic] = 
  //      cov->pref[TBM] = cov->pref[SpectralTBM] = 0;
  if ( (PisNULL(DANISO) || isMiso(mtype)) && PisNULL(DPROJ)) {
    cov->logspeed = sub->logspeed * P0(DVAR);
  }
  //////////////// 

  if (sub->pref[Nugget] > PREF_NONE) cov->pref[Specific] = 100;
  if (PisNULL(DPROJ)) cov->matrix_indep_of_x = sub->matrix_indep_of_x;

  cov->Sdollar->simplevar = simplevar;
  ASSERT_ONESYSTEM;
  cov->Sdollar->orig_owniso = owniso; // for struct Sproc !

  if (isnowProcess(cov)) {
    MEMCOPY(cov->pref, PREF_NOTHING, sizeof(pref_shorttype)); 
  } else {
    model *calling = cov->calling; 
    if (calling != NULL && isnowProcess(calling) &&
	!PisNULL(DPROJ)) {
      printf("unkom,ment spaeter"); //erst sequenial "reparieren"
      // BUG;
      //for (int i=0; i<Forbidden; i++) cov->pref[i] *= 0.5;
      //cov->pref[Specific] = 5;
    }
  }

  if (global->coords.coord_system == earth &&
      (PisNULL(DPROJ) ||
       (P0INT(DPROJ) != PROJ_TIME && (Nproj != 1 || PPROJ[0] <= 1))) && // OK
      isCartesian(DEFSYS(next)) &&//is_all(isCartesian, DefList + NEXTNR) &&
      global->messages.warn_scale &&
      (PisNULL(DSCALE) || 
       P0(DSCALE) < (STRCMP(global->coords.newunits[0], "km")== 0 ? 10 : 6.3))
      && !parallel()) {
    global->messages.warn_scale = false; // OK
    PRINTF("Value of the scale parameter equals '%4.2f', which is less than 100,\nalthough models defined on R^3 are used in the context of earth coordinates,\nwhere larger scales are expected. (This message appears only ones per session.)\n", PisNULL(DSCALE) ? 1.0 : P0(DSCALE)); 
  }
 
 // printf("checkS! end!\n");

  if ((err = TaylorS(cov)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;
}


bool allowedDS(model *cov) {
  model 
    *Aniso =cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO], // nur cartesian, nur kernel (ausser angle)
    *Daniso = cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE], // && !isRandom(Scale) if (!isDollarProc(cov) && !isKernel(OWN)) RETURN_ERR(ERRORFAILED); NUR cartesian
    *Var = cov->kappasub[DVAR];
  bool angle = isAngle(Aniso) || isAngle(Daniso),
    *D = cov->allowedD;

  if ((Scale != NULL && !isRandom(Scale) && !isDollarProc(cov) ) ||
      (Aniso != NULL && !angle) || // angle erlaubt XONLY!
      (Var != NULL && ! isRandom(Var))) {
    assert(LAST_DOMAINUSER == 1);
    D[XONLY] = false;
    D[KERNEL] = true;
    return false;
  }
  return allowedDstandard(cov);
}


bool allowedIS(model *cov) {
  model 
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],  // nur cartesian, nur kernel (auch angle)
    *Daniso = cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE], // && !isRandom(Scale) if (!isDollarProc(cov) && !isKernel(OWN)) RETURN_ERR(ERRORFAILED); NUR cartesian
    *Var = cov->kappasub[DVAR];
  bool ScaleAniso = (Scale != NULL && !isRandom(Scale)) || Aniso != NULL ||
                  Daniso != NULL,
    //angle = isAngle(Aniso),
    var = Var != NULL && !isRandom(Var);
  bool *I = cov->allowedI;

  bool allowed = allowedIstandard(cov);
  if (allowed) for (int i=FIRST_ISOUSER; i<=LAST_ISOUSER; I[i++] = true);
  //
  //for (int i=0; i<LAST_ISOUSER; i++) printf("%d ", I[i]); BUG;
  
  //  assert(!angle);
  //  printf("%d %d %d\n", Scale != NULL && !isRandom(Scale), Aniso != NULL || Daniso != NULL, !angle);
  //  APMI(cov);

  int nproj = cov->nrow[DPROJ], // OK
    *p = PINT(DPROJ);// OK
  if (nproj>0) {
    allowed = false;
    assert(SYMMETRIC == 2 + DOUBLEISOTROPIC);
    
    bool hasprev = PREV_INITIALISED;
    isotropy_type previso = hasprev ? PREVISO(0) : ISO_MISMATCH;
    bool anyEarth = hasprev && (isEarth(previso) || isSpherical(previso)),
      time = p[0] == PROJ_TIME  || (nproj == 1 && p[0] == OWNLOGDIM(0));
    if (!time && anyEarth) {      
      int min = p[0];  
      for (int i=1; i<nproj; i++) if (p[i] < min) min = p[i];
      time = min > 2;
    }
    if (time) {
      cov->IallowedDone = false;
      if (anyEarth && !isAnyIsotropic(previso))
      	set_iso(PREV, 0,
		equalsAnySymmetric(previso) ? SYMMETRIC : CARTESIAN_COORD);
      if (hasprev && (!anyEarth || PREVLOGDIM(0)<=2)) {
	for (int i = 1 + (int) CARTESIAN_COORD; i <= (int) LAST_ISOUSER;
	     I[i++]=false);
      }
    }
  }

  // aniso matrix

  if (ScaleAniso || var) {
     allowed = false;
   int i = (int) FIRST_ISOUSER,
      last = //angle ? SYMMETRIC :
      CARTESIAN_COORD;
    
    while (i <= last && !I[i]) i++;
    I[last] = i <= last;
    for ( ; i < (int) last; I[i++] = false);
    if (ScaleAniso) {
      for (i = 1 + (int)CARTESIAN_COORD; i <= (int) LAST_ISOUSER; I[i++]=false);
      return false; 
    } else {
      i = (int) FIRST_EARTH;
      while (i <= (int) EARTH_COORD && !I[i]) i++;
      I[EARTH_COORD] = i <= EARTH_COORD;
      for ( ; i < (int) EARTH_COORD; I[i++] = false);

      i = (int) FIRST_SPHERICAL;
      while (i <= (int) SPHERICAL_COORD && !I[i]) i++;
      I[SPHERICAL_COORD] = i <= SPHERICAL_COORD;
      for ( ; i < (int) SPHERICAL_COORD; I[i++] = false);
    }
  }

  else if ((!PisNULL(DANISO) /* internal */ && cov->nrow[DANISO] > 1) ||
	   (!PisNULL(DAUSER) && cov->ncol[DAUSER] > 1) || nproj > 0){
    allowed = false;
   int i = (int) FIRST_ISOUSER;
    while (i <= (int) SYMMETRIC && !I[i]) i++;
    I[SYMMETRIC] = i <= SYMMETRIC;
    for ( ; i < (int) SYMMETRIC; I[i++] = false);

    I[EARTH_SYMMETRIC] |= I[EARTH_ISOTROPIC];
    I[EARTH_ISOTROPIC] = false;
    I[SPHERICAL_SYMMETRIC] |= I[SPHERICAL_ISOTROPIC];		       
    I[SPHERICAL_ISOTROPIC] = false;

    if (LocPrevTime(cov)) {      
      I[DOUBLEISOTROPIC] = (!PisNULL(DANISO) && cov->nrow[DANISO] == 2) ||
	(!PisNULL(DAUSER) && cov->ncol[DAUSER] == 2) ||
	nproj > 0; // nicht sehr praezise...
    }   
  }

  return allowed;
}


void rangeS(model *cov, range_type* range){
  int i;
  bool negdef = isnowNegDef(cov);
  range->min[DVAR] = negdef ? 0.0 : RF_NEGINF;
  range->max[DVAR] = RF_INF;
  range->pmin[DVAR] = negdef ? 0.0 : -10000;
  range->pmax[DVAR] = 100000;
  range->openmin[DVAR] = !negdef;
  range->openmax[DVAR] = true;

  range->min[DSCALE] = 0.0;
  range->max[DSCALE] = RF_INF;
  range->pmin[DSCALE] = 0.0001;
  range->pmax[DSCALE] = 10000;
  range->openmin[DSCALE] = true;
  range->openmax[DSCALE] = true;

  for (i=DANISO; i<= DAUSER; i++) {
    range->min[i] = RF_NEGINF;
    range->max[i] = RF_INF;
    range->pmin[i] = - 10000;
    range->pmax[i] = 10000;
    range->openmin[i] = true;
    range->openmax[i] = true;
  }
  
  range->min[DPROJ] = -2;
  range->max[DPROJ] = GATTERTOTALXDIM;
  range->pmin[DPROJ] = 1;
  range->pmax[DPROJ] =  range->max[DPROJ];
  range->openmin[DPROJ] = false;
  range->openmax[DPROJ] = false;
}


Types TypeS(Types required, model *cov, isotropy_type required_iso){
  bool ok =(COVNR != DOLLAR_PROC &&
	    (isShape(required) || isTrend(required) || equalsRandom(required)))
    || (COVNR == DOLLAR_PROC && isProcess(required));
  //  printf("ok=%d\n", ok);
  if (!ok) return BadType;
 
  model *sub = cov->key==NULL ? cov->sub[0] : cov->key;

  return TypeConsistency(required, sub, required_iso);
}


void spectralS(model *cov, gen_storage *s, double *e) {
  model *next = cov->sub[DOLLAR_SUB];
  int d,
    ncol = PisNULL(DANISO) ? OWNLOGDIM(0) : cov->ncol[DANISO];
  double sube[MAXTBMSPDIM],
    *scale =P(DSCALE),
    invscale = 1.0;

  SPECTRAL(next, s, sube); // nicht gatternr

  // Reihenfolge nachfolgend extrem wichtig, da invscale auch bei aniso
  // verwendet wird


  if (scale != NULL) invscale /= scale[0];
  
  if (!PisNULL(DANISO)) {
    int 
      nrow = cov->nrow[DANISO];
    Long j, k, m,
      total = ncol * nrow;
    double
      *A = P(DANISO); 
    for (d=0, k=0; d<nrow; d++, k++) {
      e[d] = 0.0;
      for (m=0, j=k; j<total; j+=nrow) {
	e[d] += sube[m++] * A[j] * invscale;
      }
    }
  } else { 
    for (d=0; d<ncol; d++) e[d] = sube[d] * invscale;
  }

}


void ScaleDollarToLoc(model *to, model *from,
		      int VARIABLE_IS_NOT_USED depth) {

  assert(!PARAMisNULL(to, LOC_SCALE));
  assert(isDollar(from));
  assert(!PARAMisNULL(from, DSCALE));
  PARAM(to, LOC_SCALE)[0] = PARAM0(from, DSCALE);
}

bool ScaleOnly(model *cov) {
  return isDollar(cov) && 
    PisNULL(DPROJ) && cov->kappasub[DPROJ] == NULL &&
    PisNULL(DAUSER) &&  cov->kappasub[DAUSER] == NULL &&
    PisNULL(DANISO) &&  cov->kappasub[DANISO] == NULL &&
    (PisNULL(DVAR) || P0(DVAR)==1.0) && cov->kappasub[DVAR] == NULL;
}



int addScales(model **newmodel, model *calling, model *Scale, double scale) {  
  if (scale != 1.0) { 
    addModel(newmodel, LOC, calling);
    kdefault(*newmodel, LOC_SCALE, scale);
  }
  if (Scale != NULL) { // could happen that scale != 1.0 && Scale != NULL
    //                  since scale includes also scale from aniso
    if (!isnowRandom(Scale)) RETURN_ERR_COV(Scale, XERRORNONSTATSCALE);
    addModel(newmodel, LOC, calling);
    addSetDistr(newmodel, Scale->calling, ScaleDollarToLoc, true, MAXINT);
  }
  return NOERROR;
}

int addPGSLocal(model **Key, // struct of shape
		  model *shape,// shape itself
		  model *local_pts,
		  int dim, int vdim, Types frame) {
  globalparam *global = &((*Key)->base->global);
  // SEE ALSO addPGS in extremes.cc
  /// random coin und smith

  // versucht die automatische Anpassung einer PointShapeType-Funktion;
  // derzeit wird 
  // * PGS und
  // * STANDARD_SHAPE (weiss nicht mehr wieso -> coins?)
  // versucht; inklusive Init
  //
  // ruft t.w.  FillInPts auf

#define specific (nPOISSON_SCATTER - 1)
  bool maxstable = hasMaxStableFrame(shape);
  int i,
    err = NOERROR,
    method = global->extreme.scatter_method,
    pgs[specific] = {maxstable ? ZHOU : BALLANI, STANDARD_SHAPE}; 
  #define msgN 200
  char msg[nPOISSON_SCATTER - 1][LENERRMSG];
  assert(shape != NULL);
  assert(*Key == NULL);

  // to do: strokorbball: raeumlich waere Standard/Uniform besser;
  // praeferenzen programmieren?
  for (i=0; i<specific; i++) {
    if (method != POISSON_SCATTER_ANY && method != i) continue;
    
     if (i > 0) {
      errorMSG(err,shape->base, msg[i-1]);
      //     XERR(err);  // eigentlich muss das hier weg
    }
    //    if (i > 0) XERR(err); assert(i ==0);
     if (*Key != NULL) COV_DELETE(Key, shape);
    assert(shape->calling != NULL);
    addModel(Key, pgs[i], shape->calling);
 
    if ((err = FillInPts(*Key, shape)) != NOERROR) continue;
    if (MODELNR(*Key) != ZHOU) continue;
    model *local, *dummy,
      *cov = *Key;
    if ((err = covcpy(&local, false, local_pts, shape->prevloc, NULL, 
		      true, true, false)) != NOERROR) RETURN_ERR(err);
    SET_CALLING(local, (*Key)->calling);
    dummy = local;
    while (dummy->sub[DOLLAR_SUB] != NULL) dummy = dummy->sub[DOLLAR_SUB];
    if (MODELNR(dummy) != LOC) BUG;
    dummy->sub[DOLLAR_SUB] = *Key;
    SET_CALLING(*Key, dummy);
    
    SET_CALLING(cov, shape->calling);
    SET_CALLING(cov->sub[PGS_FCT], cov);
    SET_CALLING(cov->sub[PGS_LOC], cov);

    assert(cov->sub[PGS_LOC] != NULL && cov->sub[PGS_FCT] != NULL);
    cov->nsub = 2;

     if ((err = CHECK(*Key, dim, dim, PointShapeType, XONLY, 
		     CoordinateSystemOf(ISO(SYSOF(shape), 0)),
		     vdim, frame)) != NOERROR) {
      continue; 
    }
    ONCE_NEW_COV_STORAGE(cov, gen);

    if ((err = INIT(cov, 1, cov->Sgen)) == NOERROR) break;
  } // for i_pgs

  model *cov = *Key;
  if (err != NOERROR) {
     SERR("error occured when creating the local point-shape fctn");
   // errorstring_type xx = "error occured when creating a local point-shape fctn";
    // errorMSG(err, xx, cov->base, msg[i-1]);
    // SERR2("no method found to simulate the given model:\n\tpgs: %.50s\n\tstandard: %.50s)", msg[0], msg[1]);
  }
  RETURN_NOERROR;
}



int structS(model *cov, model **newmodel) {
  if (isnowProcess(cov)) {
    assert(hasGaussMethodFrame(cov));
    SET_NR(cov, DOLLAR_PROC);
    return structSproc(cov, newmodel); // kein S-TRUCT(...) !!
  } 

  model
    *next = cov->sub[DOLLAR_SUB],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *Scale = cov->kappasub[DSCALE];
  int err = NOERROR; 
  bool generalAniso = (Aniso != NULL && !isAngle(Aniso)) ||
    (Scale != NULL && !isnowRandom(Scale));

  if (generalAniso)
    SERR1("complicated models including arbitrary functions for '%.50s' cannot be simulated yet", KNAME(DAUSER));
   
  ASSERT_NEWMODEL_NOT_NULL;
 
  if (cov->kappasub[DVAR] != NULL) {
    SERR2("Arbitrary functions for '%.50s' should be replaced by multiplicative models using '%.50s'", KNAME(DVAR), DefList[PROD].nick);
  } 

  switch (cov->frame) {
    // case Poisson : find location scattering to given shape
  case SmithType : { // distribution for locations for smith
    // only for compactly supported shape functions
    ASSERT_ONESYSTEM; 
    ASSERT_NEWMODEL_NOT_NULL; // ?
    if (next->randomkappa) SERR("arbitrary random shapes not programmed yet");
    if (!PisNULL(DANISO))RETURN_ERR(ERRORANISO);
    if (!PisNULL(DANISO)) {
      if (!isMonotone(next) || !isIsotropicXonly(next) ||
	  PisNULL(DANISO) || cov->ncol[DANISO] != cov->nrow[DANISO])
	SERR("anisotropy parameter only allowed for simple models up to now");
    }
    
    if (!PisNULL(DPROJ)) SERR("projections not allowed yet");
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);
    if (*newmodel == NULL) // i.e. no information about how to model
      SERR("probably not a compactly supported shape function");
    
    assert(*newmodel != NULL);
   
    double anisoScale = 1.0;  
    if (!PisNULL(DANISO)) {
      anisoScale = 1.0 / getMinimalAbsEigenValue(P(DANISO), cov->nrow[DANISO]);
      if (!R_FINITE(anisoScale)) RETURN_ERR(ERRORANISO_INV);
      if (anisoScale < 0) RETURN_ERR(-anisoScale);
    }
    
      
    Types type = BadType;
    double scale = PisNULL(DSCALE) ? 1.0 :P0(DSCALE);
  
    int
      xdim = OWNXDIM(0),
      logdim = OWNLOGDIM(0);
     
    // in same cases $ must be set in front of a pointshape function
    // and becomes a pointshape function too, see addpointshapelocal here
    if (isPointShape(*newmodel)) type = PointShapeType;
    else if ((err = CHECK_R(*newmodel, logdim)) == NOERROR) {
      type = RandomType;
      if ((err=addScales(newmodel, *newmodel, Scale, scale * anisoScale))
	  !=NOERROR) RETURN_ERR(err);
      SET_CALLING(*newmodel, cov);      
    } 

    if (!equalsRandom(type)) {
      type = type == BadType ? ShapeType : type;
      if ((err = CHECK(*newmodel, logdim, xdim, type, OWNDOM(0), OWNISO(0),
		       cov->vdim, cov->frame)) != NOERROR) RETURN_ERR(err);

      // ??   addModel(newmodel, FIRSTDOLLAR); // 2.2.19
      // ??assert( (*newmodel)->calling == cov);

      // ??  if (!PisNULL(DVAR)) kdefault(*newmodel, DVAR, P0(DVAR));

      // ?? double scale = 1.0;
      //    if (!PisNULL(DSCALE)) {
      //      kdefault(*newmodel, DSCALE, P0(DSCALE));
      //      scale = P0(DSCALE);}
 
      if (type == PointShapeType && 
	  (err = addScales((*newmodel)->sub + PGS_LOC, *newmodel, Scale,
			   anisoScale * scale)) != NOERROR) RETURN_ERR(err);
      if ((err = CHECK(*newmodel, logdim, xdim, type, PREVDOM(0),
		       PREVISO(0), cov->vdim, cov->frame)) != NOERROR)
	RETURN_ERR(err);
      if (equalsPointShape(type)) {
	BUG;
	// i.e. the case where scale calls Huetchen.cc
	// when does this happen????
	if ((err = FillInPts(*newmodel, cov)) != NOERROR) RETURN_ERR(err);
      } else {
	model *local = NULL,
	  *dummy = *newmodel;
	assert(dummy != NULL);
	*newmodel = NULL;
	// suche nach geeigneten locationen
	if ((err = addScales(&local, dummy, Scale, scale * anisoScale))
	    != NOERROR) {
	  COV_DELETE(&dummy, cov);
	  if (local != NULL) COV_DELETE(&local, cov);
	  RETURN_ERR(err)
	}
	if ((err = addPGSLocal(newmodel, dummy, local, OWNLOGDIM(0), VDIM0,
			       cov->frame)) != NOERROR) {
	  COV_DELETE(&dummy, cov);
 	  assert(local != NULL);
	  COV_DELETE(&local, cov);
	  RETURN_ERR(err);
	}
      }
    } // ! randomtype
    }
  break;
  
  case SchlatherType: case BrMethodType:
    if (next->randomkappa) SERR("random shapes not programmed yet");
    if (!PisNULL(DPROJ)) SERR("projections not allowed yet");
    // P(DVAR) hat keine Auswirkungen
    if (!PisNULL(DANISO)) {
      if (!isMonotone(next) || !isIsotropicXonly(next) ||
	  PisNULL(DANISO) || cov->ncol[DANISO] != cov->nrow[DANISO])
	SERR("anisotropy parameter only allowed for simple models up to now");
    }
    
    assert(cov->calling != NULL); 
    
    // BUG;
    // 
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);   
    break;
    
  case GaussMethodType :
  case PoissonGaussType: //find shape to given covariance function
    if (cov->key != NULL) COV_DELETE(&(cov->key), cov);//unnoetig bei PoissonG 

    if (LocDist(cov)) 
      SERR("distances do not allow for more sophisticated simulation methods");
    
    if ((err = STRUCT(next, newmodel)) > NOERROR) RETURN_ERR(err);

    assert(*newmodel != NULL);

    //    PMI(*newmodel); PMI(cov);

    if ( // next->randomkappa || -- leider erst ab check bekannt!
	!PisNULL(DPROJ) || !PisNULL(DANISO) ||
	!PisNULL(DAUSER) || cov->kappasub[DANISO] != NULL ||
	cov->kappasub[DSCALE] != NULL || cov->kappasub[DVAR] != NULL
	 )
      SERR("random shapes and more complicated shape functions not programmed yet");

    if (isRandom(*newmodel)) {
      addModelX(newmodel, LOC);
      if (!PisNULL(DVAR) && P0(DVAR) != 1.0)
	ERR("variance of the random location cannot be changed");
      if (!PisNULL(DSCALE) ) kdefault(*newmodel, LOC_SCALE, 1.0 /
				      P0(DSCALE));
    } else {      
      addModelX(newmodel, FIRSTDOLLAR);
      if (!PisNULL(DVAR)) {
	if (hasPoissonGaussFrame(cov)) kdefault(*newmodel, DVAR, P0(DVAR));
	else kdefault(*newmodel, DVAR, SQRT(P0(DVAR)));      
      }
      if (!PisNULL(DSCALE) ) kdefault(*newmodel, DSCALE, P0(DSCALE));
    }
    
    break;
 
  default :
    BUG;
    SERR2("%.50s : changes in scale/variance not programmed yet for '%.50s'", 
	  NICK(cov), TYPE_NAMES[cov->frame]);      
  }
    
  RETURN_NOERROR;
}




int initS(model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!

   model *next = cov->sub[DOLLAR_SUB],
    *Var = cov->kappasub[DVAR],
    *Scale = cov->kappasub[DSCALE],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO],
    *projM = cov->kappasub[DPROJ];
  int 
    vdim = VDIM0,
    nm = cov->mpp.moments,
    nmvdim = (nm + 1) * vdim,
    err = NOERROR;
  bool 
    angle = isAngle(Aniso),
    smith = hasSmithFrame(cov);// Realisationsweise 

  if (smith || hasAnyPoissonFrame(cov)) {//Average !!ohne smith selbst!!
    ASSERT_ONESYSTEM;
    double
      var[MAXMPPVDIM],
      scale = PisNULL(DSCALE)  ? 1.0 : P0(DSCALE);
    int dim = OWNLOGDIM(0);
 
    if (!PisNULL(DPROJ) || !PisNULL(DANISO) || 
	projM!= NULL || (Aniso != NULL && (!angle || Aniso->randomkappa))){

      SERR("(complicated) anisotropy ond projection not allowed yet in Poisson related models");
    }
  
     // Achtung I-NIT_RANDOM ueberschreibt mpp.* !!
    if (Var != NULL) {
      int nm_neu = nm == 0 && !smith ? 1 : nm;
      if ((err = INIT_RANDOM(Var, nm_neu, s, P(DVAR))) != NOERROR) RETURN_ERR(err); 
      
      int nmP1 = Var->mpp.moments + 1;
      for (int i=0; i<vdim; i++) {
	int idx = i * nmP1;
	var[i] = smith ? P0(DVAR) : Var->mpp.mM[idx + 1];      
      }
    } else for (int i=0; i<vdim; var[i++] = P0(DVAR));

    if (Scale != NULL) {
      if (dim + nm < 1) SERR("found dimension <= 0");
      int dim_neu = smith ? nm : (dim + nm) < 1 ? 1 : dim + nm; 
      if ((err = INIT_RANDOM(Scale, dim_neu, s, P(DSCALE)))
	  != NOERROR) RETURN_ERR(err);
      scale = smith ? P0(DSCALE) : Scale->mpp.mM[1];      
    }
    if ((err = INIT(next, nm, s)) != NOERROR) RETURN_ERR(err);


    for (int i=0; i < nmvdim; i++) {
      cov->mpp.mM[i] = next->mpp.mM[i]; 
      cov->mpp.mMplus[i] = next->mpp.mMplus[i]; 
    }

    if (Var != NULL && !smith) {
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= Var->mpp.mM[i];
	cov->mpp.mMplus[i] *= Var->mpp.mMplus[i];
      }
    } else {
      int j, k;
      double pow_var;
      for (k=j=0; j<vdim; j++) { 
	pow_var = 1.0;
	for (int i=0; i <= nm; i++, pow_var *= var[j], k++) {	
	  cov->mpp.mM[k] *= pow_var;
	  cov->mpp.mMplus[k] *= pow_var;
	}
      }
    }
    if (Scale != NULL && !smith) {
      if (Scale->mpp.moments < dim) SERR("moments can not be calculated.");
      int j, k,
	nmP1 = Scale->mpp.moments + 1;
      for (k=j=0; j<vdim; j++) { 
	double pow_scale = Scale->mpp.mM[dim + j * nmP1];
	for (int i=0; i <= nm; i++, k++) {
	  cov->mpp.mM[k] *= pow_scale;
	  cov->mpp.mMplus[k] *= pow_scale;
	}
      }
    } else {
      double pow_scale = POW(scale, dim);
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= pow_scale;
	cov->mpp.mMplus[i] *= pow_scale;
      }
    }

    if (!PisNULL(DANISO)) {      
      if (cov->nrow[DANISO] != cov->ncol[DANISO]) RETURN_ERR(ERRORANISO_SQUARE);
      double invdet = FABS(1.0 / getDet(P(DANISO), cov->nrow[DANISO]));
      if (!R_FINITE(invdet)) RETURN_ERR(ERRORANISO_DET);
      
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= invdet;
	cov->mpp.mMplus[i] *= invdet;
      }

    }

    if (Aniso != NULL) {  
      double invdet;
      if (angle) {
	int 
	  ncol = Aniso->vdim[0],
	  nrow = Aniso->vdim[1];
	if (nrow != ncol) RETURN_ERR(ERRORANISO_SQUARE);
	double 
	  *diag = PARAM(Aniso, ANGLE_DIAG);
	if (diag != NULL) {
	  invdet = 1.0;
	  for (int i=0; i<ncol; i++) invdet /= diag[i];
	} else {
	  invdet = PARAM0(Aniso, ANGLE_RATIO);
	}
      } else {
	SERR("only anisotropy matrices basesd on RMangle allowed.");
      }
      
      invdet = FABS(invdet);      
      if (!R_FINITE(invdet)) RETURN_ERR(ERRORANISO_DET);
      
      for (int i=0; i < nmvdim; i++) {
	cov->mpp.mM[i] *= invdet;
	cov->mpp.mMplus[i] *= invdet;
      }
    }

    bool unnormed = R_FINITE(next->mpp.unnormedmass);
    if (unnormed) {
      //printf("unnormedmass in $ %.50s\n", NAME(next));
      if (vdim > 1) BUG;
      cov->mpp.unnormedmass = next->mpp.unnormedmass * var[0];
    } else cov->mpp.unnormedmass = RF_NAN;
    
    if (isnowRandom(next)) {
      if (unnormed) {
	int maxv = MIN(vdim, MAXMPPVDIM);
	for (int i=0; i<maxv; i++)  cov->mpp.maxheights[i] = RF_INF;
      } else RETURN_ERR(ERRORNOTPROGRAMMEDYET);
    } else {
      int maxv = MIN(vdim, MAXMPPVDIM);
      for (int i=0; i<maxv; i++) 
	cov->mpp.maxheights[i] = next->mpp.maxheights[i] * var[i]; // maxv
    }
  } // hasSmithFrame

  else if (hasGaussMethodFrame(cov) || hasAnyEvaluationFrame(cov)) {
    model 
      *key = cov->key,
      *sub = key == NULL ? next : key;
    assert(sub != NULL);
    assert(key == NULL || ({PMI(cov);false;}));//
   
    if ((err=INIT(sub, 0, s)) != NOERROR) RETURN_ERR(err);
       
  } else { 
    if ((err=INIT(next, 0, s)) != NOERROR) RETURN_ERR(err); // e.g. from MLE
  }

  if ((err = TaylorS(cov)) != NOERROR) RETURN_ERR(err);

  RETURN_NOERROR;

}


void doS(model *cov, gen_storage *s){

  //  printf("hier doS\n");
  
  model
    *Var = cov->kappasub[DVAR],
    *Scale = cov->kappasub[DSCALE];
  int vdim = VDIM0;
   
  if (Var != NULL) {
    if (isnowRandom(Var)) {
      if (isProcess(Var)) BUG;
      assert(!PisNULL(DVAR));
      DORANDOM(Var, P(DVAR));
    } else if (Var->randomkappa) {
      assert(!PisNULL(DVAR));
      DO(Var, s);
    } else BUG;
  }
  
 if (Scale != NULL) {
   if (isnowRandom(Scale)) {
      if (isProcess(Scale)) BUG;
      assert(!PisNULL(DSCALE));
      DORANDOM(Scale, P(DSCALE));
   } else if (Scale->randomkappa) {
      BUG;
      assert(!PisNULL(DSCALE));
      DO(Scale, s);
    } else BUG;
  }

  
 if (hasSmithFrame(cov) || hasAnyPoissonFrame(cov) ||
     hasEvaluationFrame(cov)) // truncsupport creates evaluation frame
   {
    model *next = cov->sub[DOLLAR_SUB];
    
    DO(next, s);// nicht gatternr
    int maxv = MIN(vdim, MAXMPPVDIM);
    for (int i=0; i<maxv; i++) 
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * P0(DVAR); // maxv
    return;
  }
  
  else if (hasGaussMethodFrame(cov)) {    
    double 
      *res = cov->rf,
      sd = SQRT(P0(DVAR));
    int 
      totalpoints = Loctotalpoints(cov);
    assert(res != NULL);
    if (cov->key == NULL) BUG;

    if (sd != 1.0) for (int i=0; i<totalpoints; i++) res[i] *= (double) sd;
    
    return;
  }

  BUG;
}

//////////////////////////////////////////////////////////////////////
// PROCESSES
//////////////////////////////////////////////////////////////////////

int structSproc(model *cov, model **newmodel) {
  
  model
    *next = cov->sub[DOLLAR_SUB],
    *Scale = cov->kappasub[DSCALE],
    *Aniso = cov->kappasub[DANISO] == NULL ? cov->kappasub[DAUSER]
    : cov->kappasub[DANISO];
  int
    tsdim = Loctsdim(cov), 
    newtsdim = tsdim,
    err = NOERROR; 
  // model *sub;
  assert(isDollarProc(cov));

  assert(newtsdim  > 0);
  
  if ((Aniso != NULL && Aniso->randomkappa) ||
      (Scale != NULL && Scale->randomkappa)
      ) {
    SERR1("complicated models including arbitrary functions for '%.50s' cannot be simulated yet", KNAME(DAUSER));
  }
  switch (cov->frame) {
  case GaussMethodType : {

    //    PMI0(cov);
    
    ASSERT_NEWMODEL_NULL;
    if (cov->key != NULL) COV_DELETE(&(cov->key), cov);

    if (LocDist(cov)) 
      SERR("distances do not allow for more sophisticated simulation methods");

    bool grid = Locgrid(cov),
      Time = LocTime(cov);
    int spatialpoints = Locspatialpoints(cov);
    coord_type gr = Locxgr(cov);
    cov->Sdollar->somegrid = cov->Sdollar->blockmatrix =
      cov->Sdollar->separable = false;
    if (Aniso!= NULL) { // Aniso fctn
      // printf("dim = %d\n", tsdim);
      assert(tsdim > 0);
      if (cov->ownloc != NULL) BUG;
       double *xx=NULL;
      Long total = Loctotalpoints(cov);
      int xdim = TransformLoc(cov, Loc(cov), &xx, NULL, True);
      //printf("back hier\n");
    
      newtsdim = Aniso->vdim[0];
      assert(LocSets(cov) == 1);
      if (empty_loc_set(cov, newtsdim, total, 0)){
	FREE(xx);
	ERR("memory allocation error when expanding Aniso definition");
      }
      int bytes = newtsdim * sizeof(double);
      double
	*x = xx, 
	*xnew = Locx(cov);
      
      DEFAULT_INFO(info);
      for (Long i=0; i<total; i++, x+=xdim, xnew += newtsdim) {
	info[INFO_IDX_X] = i;	
	FCTN(x, info, Aniso, xnew);
      }
      FREE(xx);
    } else if (Scale != NULL) {
      // check code wether it will still work
      if (isnowRandom(Scale)) {
	SERR("no specific simulation technique available for random scale");
      }
      SERR("no specific simulation technique available for arbitrary scale");
    } else {

       double
	*aniso = P(DANISO),
	scale = PisNULL(DSCALE) ? 1.0 : P0(DSCALE);
      int *proj =  PPROJ,
	nproj = Nproj,
	ncol = cov->ncol[DANISO],
	nrow = cov->nrow[DANISO];
      
      cov->Sdollar->separable =
	Loc(cov)->caniso == NULL && (grid || Time);

      // printf("Pisn %d \n", PisNULL(DPROJ));

      // && isXonly(NEXT); 21.1.21: diese Bedingungen entfernt,
      //                            entsprechend gewuenscht(!!) eingeschraenkt
      //                            sein soll Siehe Bsp RMfix.Rd
      // cov model may (!!) depend on the location!
      //                insb. wenn i_row/i_col verwendet wurde, muss
      //                kernel irgendwo dazwischen auftauchen
      cov->Sdollar->separable = !PisNULL(DPROJ);
      if (!cov->Sdollar->separable &&
	  !PisNULL(DANISO) && cov->ncol[DANISO] < cov->nrow[DANISO]
	  // TO DO dies kann auch nochmals abgeschwaecht werden
	  // dann muss aber sehr auf raum-zeit aufgepasst werden
	  ) {
	// possibly 'aniso' is a projection matrix 
	// in the widest sence (instead of +1 any value != 0)	
	bool truelySpatial = !grid && ncol > 1;

	if ((cov->Sdollar->separable = isMproj(Type(aniso, nrow, ncol)))) {
	  if (truelySpatial) {
	    for (int d=0; d<ncol; d++)
	      if (aniso[d * nrow - 1] != 0.0) { // i.e. there is time comp.
		cov->Sdollar->separable = false;
		break;
	      }
	  }
	  
	  if (cov->Sdollar->separable) { // separability still not sure
	    // projections must be in different directions
	    nproj = Nproj = ncol;
	    int bytes = sizeof(int) * nproj;
	    proj = PPROJ = (int*) MALLOC(bytes);
	    bool *used = (bool*) CALLOC(nproj, sizeof(bool));	 
	    for (int i=0; i<ncol; i++) {
	      double *a=aniso + i * nrow;
	      for (int j=0; j<nrow; j++) {
		if (a[j] != 0.0) {
		  cov->Sdollar->separable &= !used[j];		  
		  used[j] = true;
		  proj[i] = j;
		  break;
		  }
	      }
	      if (!cov->Sdollar->separable) break;
	    }
	    FREE(used);
	  }
	}
      } // !PisNULL(aniso)
 
      if (cov->Sdollar->separable) {  	
	ALLC_NEWINT(Sdollar, len, tsdim  + 1, len);
	ALLC_NEWINT(Sdollar, cumsum, tsdim + 1, cumsum);
	ALLC_NEWINT(Sdollar, total, tsdim + 1, total);
	
	double *T= LocT(cov);   
	for (int d=0; d<tsdim; d++) {
	  total[d] = cumsum[d] = 0;
	  if (grid) len[d] = (int) gr[d][XLENGTH];
	  else if (d == 0) len[d] = spatialpoints;
	  else if (d == tsdim - 1 && Time) len[d] = (int) T[XLENGTH];
	  else len[d] = 1;
	}
	
	cov->Sdollar->somegrid =
	  grid || (nproj == 1 && Time && proj[0] == tsdim - 1);
	if (cov->Sdollar->somegrid) {
	  double *xx = (double*) MALLOC(sizeof(double) * tsdim * 3);
	  for (int d=0; d<nproj; d++) {
	    int neu = proj[d],
	      alt = proj[d-1];
	    MEMCOPY(xx + d * 3, neu == tsdim - 1 && Time ? T : gr[neu],
		    3 * sizeof(double));
	    if (aniso == NULL) xx[d * 3 + XSTEP] /= scale;
	    else xx[d * 3 + XSTEP] *= aniso[d * nrow + neu];
	    if (d == 0) cumsum[neu] = 1;
	    else cumsum[neu] = cumsum[alt] * (int) gr[alt][XLENGTH];
	    total[neu] = cumsum[neu] * (int) gr[neu][XLENGTH];
	  }
	  //	  printf("some grdi\n")	 ; 
	  loc_set(xx, NULL, NULL, NULL, nproj, nproj, UNKNOWN_NUMBER_GRIDPTS, 0,
		  false, true, true, false, cov);
	  FREE(xx);
	} else { // purely space, but can can coordinates might be permuted
	  for (int i=0; i<nproj; i++) // check if purely spatial
	    cov->Sdollar->separable &= proj[0] != tsdim - 1;	      
 	  if (cov->Sdollar->separable) {
	    int xdimOZ = LocxdimOZ(cov);	    
	    Long bytes = (Long) sizeof(double) * nproj * spatialpoints;
	    double *x0 = (double*) MALLOC(bytes),
	      *xx = x0,
	      *A = x0 + nproj * (spatialpoints - 1), // einfach nur ein
	    // dummy speicherort, der exakt zum Schluss ueerschrieben wird
	      *L = Locx(cov);
	    if (aniso != NULL)
	      for (int d=0 ; d<nproj; d++) {
		A[d] = aniso[d * nrow +  proj[d]];
		assert(A[d] != 0.0);
	      }					     
	    //  printf("AFE aniso = %d %f\n", aniso!=NULL, scale);
	    for (int i=0 ; i<spatialpoints; i++, L += xdimOZ) {
	      //    printf("i=%d L=%ld\n", i, L-Locx(cov));
	      for (int d=0 ; d<nproj; d++) {
		//		if (i > 105) printf("  d=%d %d\n", d, proj[d]);
		if (aniso == NULL) *(xx++) = L[proj[d]] /= scale;
		else *(xx++) = L[proj[d]] * A[d];
	      }
	    }
	    //  	    printf("no grdi %d %d %d %d\n",nproj,spatialpoints,grid,PPROJ[0]); 
	    loc_set(x0, NULL, NULL, NULL, nproj, nproj, spatialpoints, 0,
		    false, grid, false, false, cov);
	    FREE(x0);
	  }
	}
      }
 
      if (!cov->Sdollar->separable) {
	// not a projection matrix; Block matrices with block of
	// zeros are still separable, spatially, or if it is on a grid
	cov->Sdollar->separable = false;
	// to do
	// cov->Sdollar->blockmatrix = true;
	// could stem from a projection matrix with multiple projection onto
	// the same direction
      }
      
      if (!cov->Sdollar->separable) {
	usr_bool gridexpand = NEXTNR == TBM_PROC_INTERN || isAnyNugget(NEXTNR)
	  ? False : GRIDEXPAND_AVOID;
	TransformLoc(cov, true /* timesep*/, gridexpand,  // OK
		       True /* involveddollar */);
      }
      newtsdim = Loctsdim(cov);
      assert(newtsdim > 0 );
    } // Scale and Aniso are not given as functions


    //   printf("newtsdim = %d %d\n",newtsdim, cov->Sdollar->separable);

   //PMI0(cov);
    //   TREE(cov);

    
    if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
    if (!isGaussMethod(cov->key)) addModelKey(cov, GAUSSPROC);
    SetLoc2NewLoc(cov->key, LocP(cov));
    
    model *key;
    key = cov->key;
    assert(key->calling == cov);
    assert(newtsdim > 0);
    
    ASSERT_ONESYSTEM;
    if ((err = CHECK_NO_TRAFO(key, newtsdim, newtsdim, ProcessType, XONLY, 
			      CoordinateSystemOf(cov->Sdollar->orig_owniso),
			      VDIM0, GaussMethodType)) != NOERROR) {
      //APMI(cov);
      RETURN_ERR(err);
    }

    //PMI(cov);
    
    err = STRUCT(key, NULL);

    RETURN_ERR(err);
  }
  default :
    SERR2("%.50s: changes in scale/variance not programmed yet for '%.50s'", 
	  NICK(cov), TYPE_NAMES[cov->frame]);      
  }
     
  RETURN_NOERROR;
}




int initSproc(model *cov, gen_storage *s){
  // am liebsten wuerde ich hier die Koordinaten transformieren;
  // zu grosser Nachteil ist dass GetDiameter nach trafo 
  // grid def nicht mehr ausnutzen kann -- umgehbar?!

  model *key = cov->key;
  location_type *prevloc = LocPrev(cov);
  int 
    prevdim = prevloc->timespacedim,
    prevtotalpts = prevloc->totalpoints,
    owntotalpts =  Loctotalpoints(cov);

  int err = NOERROR;

  assert(key != NULL);
  
  if ((err = INIT(key, 0, s)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  key->simu.active = true; 
  assert(s != NULL);
  cov->fieldreturn = wahr;
  cov->origrf = prevtotalpts != owntotalpts;
  assert(cov->Sextra != NULL);
  if (cov->origrf) {
    assert(prevtotalpts % owntotalpts == 0);
    assert(!PisNULL(DPROJ) || !PisNULL(DANISO));
    assert(VDIM0 == VDIM1);
    // projection or reducing anisotropy matrix
    //    printf("\n\n\n\n\n ********* cov-rf %d \n",  VDIM0 * prevtotalpts);
    cov->rf = (double*) MALLOC(sizeof(double) * VDIM0 * prevtotalpts);
  } else cov->rf = cov->key->rf;

  model *Var = cov->kappasub[DVAR];
  if (Var != NULL) {
    if (isnowRandom(Var) || Var->randomkappa) {
      if (isProcess(Var)) {
	RETURN_ERR(ERRORNOTPROGRAMMEDYET);
      } else {
	if ((err = INIT(Var, 0, s)) != NOERROR)
	  RETURN_ERR(err);	  
      }
    } else {
      assert(cov->Sdollar != NULL);
      int totptsvdim = prevtotalpts * VDIM0;
      cov->Sdollar->sd = (double*) MALLOC(sizeof(double) * totptsvdim);
      double *sd = cov->Sdollar->sd;
      FctnExtern(cov, Var, true, sd);
      for (int i=0; i<totptsvdim; i++) sd[i] = SQRT(sd[i]);
    }
  }
  
  model *Scale = cov->kappasub[DSCALE];
  if (Scale != NULL) {
    RETURN_ERR(ERRORNOTPROGRAMMEDYET);
  }
 
  RETURN_NOERROR;
}


//int zz = 0;
void doSproc(model *cov, gen_storage *s){
  int 
    vdim = VDIM0; 

  if (hasGaussMethodFrame(cov)) {    
    assert(cov->key != NULL);
    double *res = cov->key->rf;
    int 
      totalpoints = Loctotalpoints(cov),
      totptsvdim = totalpoints * vdim;

    //    for (int i=0; i<VDIM0 * LocLoctotalpoints(LocPrev(cov)); i++) cov->rf[i] = i;
    //    printf("short do %s %e %e %d Time=%d\n", NAME(cov), cov->rf[0], cov->rf[VDIM0 * LocLoctotalpoints(LocPrev(cov))-1], LocLoctotalpoints(LocPrev(cov)), LocLoctime(LocPrev(cov)));  return; 
   
    DO(cov->key, s);

    //    printf("hiere\n");   
    
    model *Var = cov->kappasub[DVAR];

    if (Var != NULL) {
      if (isnowRandom(Var) || Var->randomkappa) {
	if (isProcess(Var)) {
	  XERR(ERRORNOTPROGRAMMEDYET);
	} else {
	  DORANDOM(Var, P(DVAR));
 	}
	double sd = SQRT(P0(DVAR));
	for (int i=0; i<totptsvdim; i++) res[i] *= sd;
 	
      } else {
	double *sd = cov->Sdollar->sd;
	assert(sd != NULL);
	for (int i=0; i<totptsvdim; i++) res[i] *= sd[i];
      }
    } else {
      assert(!PisNULL(DVAR));
      double sd = SQRT(P0(DVAR)); 
      if (sd != 1.0) {
	for (int i=0; i<totptsvdim; i++) res[i] *= sd;
      }
    }
  }
  
  else if (hasMaxStableFrame(cov) || hasAnyPoissonFrame(cov)) {
    BUG; // 2.2.19: darf eigentlich nicht (so) genutzt werden
    assert(vdim == 1);
    model *next = cov->sub[DOLLAR_SUB],
      *Var = cov->kappasub[DVAR],
      *Scale = cov->kappasub[DSCALE];
 
    assert(VDIM0 == VDIM1);

    if (Var != NULL && Var->randomkappa) {
      BUG; // Da muss irgenduas dann auf rf draufmultipliziert werden
      assert(!PisNULL(DVAR) && isnowRandom(Var));
      DEFAULT_INFO(info);
      VTLG_R(NULL, info, Var, P(DVAR));
    }

    if (Scale != NULL && Scale->randomkappa) {
      // remote could be deterministic although local ist not
      BUG; // otherwise do must depend on the new scale oder man muss
      // interpolieren o.ae.
      assert(!PisNULL(DSCALE));
      DEFAULT_INFO(info);
      VTLG_R(NULL, info, Scale, P(DSCALE));
    }
    
    DO(next, s);// nicht gatternr

    int maxv = MIN(vdim, MAXMPPVDIM);
    for (int i=0; i<maxv; i++)   
      cov->mpp.maxheights[i] = next->mpp.maxheights[i] * P0(DVAR); // maxv
  }

  else {
    //PMI(cov);
    BUG;
  }
 
  if (cov->origrf) { // es wird in cov->key->rf nur 
    // der projezierte Teil simuliert. Dieser Teil
    // muss auf das gesamte cov->rf hochgezogen werden
    //if (vdim != 1) BUG;
    assert(cov->Sdollar->separable)
    location_type *prevloc = LocPrev(cov);
    int
      prevtotalpts = prevloc->totalpoints,
      owntotalpts =  Loc(cov)->totalpoints,
      owntotptsvdim =  owntotalpts * VDIM0;
     
    if (cov->Sdollar->blockmatrix) {
    }
    
    else if (cov->Sdollar->separable) {
      assert(cov->key != NULL);

      if (cov->Sdollar->somegrid) {
	int prevdim = prevloc->timespacedim;
	TALLOC_L1(nx, prevdim);      
	for (int d=0; d<prevdim; d++) nx[d] = 0;

	//	printf("%d %d %d %d\n", cov->Sdollar->somegrid, prevdim, 1, 1);
      
	
	int d=0,
	  zaehler = 0,
	  i = 0,
	  *cumsum = cov->Sdollar->cumsum,
	  *len = cov->Sdollar->len,
	  *total = cov->Sdollar->total;
	assert(total != NULL && cumsum != NULL && len != NULL);
	
	for (int v=0; v<vdim; v++) {
	  double *res = cov->rf + v * prevtotalpts,
	    *rf = cov->key->rf + v * owntotalpts;
	  while (true) {
	    res[i++] = rf[zaehler];
	    d = 0;			
	    nx[d]++;			
	    zaehler += cumsum[d];
	    while (nx[d] >= len[d]) {	
	      nx[d] = 0;		
	      zaehler -= total[d];
	      if (++d >= prevdim) break;	
	      nx[d]++;			
	      zaehler += cumsum[d];					
	    }
	    if (d >= prevdim) break;			
	  }
	}
      
	END_TALLOC_L1;

	//	BUG; printf("hier\n");
	
      } else { // purely spatial
	//PMI0(cov);
 	assert(LocLocspatialpoints(prevloc) == Loctotalpoints(cov));
	assert(LocLocTime(prevloc));
	assert(LocLocT(prevloc) != NULL)
	int
	  bytes = owntotptsvdim * sizeof(double),
	  n = (int) LocLocT(prevloc)[XLENGTH];
	double *res = cov->rf;
	assert(res != NULL);
	assert(cov->key != NULL);
	assert(cov->key->rf != NULL);

	for (int i=0; i<n; i++, res += owntotptsvdim) {
	  MEMCOPY(res, cov->key->rf, bytes);
	}
      }
    }
    
    else {
      assert(prevtotalpts == owntotalpts);
    }
   
  }
  //  printf("%d\n", LocPrev(cov)->totalpoints * vdim);
  //  printf("$proc: %e %e\n", cov->rf[0], cov->rf[LocPrev(cov)->totalpoints * vdim-1]);
}

