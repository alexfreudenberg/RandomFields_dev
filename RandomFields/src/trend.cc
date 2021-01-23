
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2011 -- 2015 Marco Oesting & Martin Schlather
               2015 -- 2017 Martin Schlather


 Handling of the different possibilities to pass the trend

Note:
 * Never use the below functions directly, but only by the functions indicated 
   in RFsimu.h, since there is no error check (e.g. initialization of RANDOM)
 * VARIANCE, SCALE are not used here 
 * definitions for the random coin method can be found in MPPFcts.cc
 * definitions for genuinely anisotropic or nondomain models are in
   SophisticatedModel.cc; hyper models also in Hypermodel.cc
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

/*Note:
the parameter 'polycoeff' is hidden in R; it is a vector consisting of the
coefficients of the corresponding trend polynomials:
the first choose(polydeg[1]+d,d) components belong to the first polynomial,
the next choose(polydeg[2]+d,d) to the second one,...
the corresponding monomial functions are of the following order:
1, x, x^2, ..., x^k, y, x y, ... x^(k-1) y, y^2, x y^2..., y^k,
z, x z, x^2 z, ...., x^(k-1) z, y z, x y z, x^2 y z, ..., z^k
*/


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

#include "questions.h"
#include "Processes.h"
#include "shape.h"
#include "Coordinate_systems.h"
#include "rf_interfaces.h"
#include "QMath.h"



int binomialcoeff(int n, int k) {
  //programmed as in wikipedia
  int i, res=1;
  
  if((k < 0) || (k > n)) return 0;
  if(k > n-k) k = n-k; //symmetry
  for(i=0; i<k; i++) {
     res *= n-i;
     res /= i+1; //two steps because of integer division
  }
  return res;
}



//////////////////////////////////////////////////////////////////////
//    shape
//////////////////////////////////////////////////////////////////////

void shapefct(double *x, int *info, model *cov, double *v){
  int vdim = VDIM0;

  assert(isnowShape(cov) || isnowTrend(cov));
  if (hasAnyEvaluationFrame(cov)) {
    int vSq = vdim * vdim;
    //    PMI0(cov);
    // BUG;
    for (int i=0; i<vSq; i++) v[i]=0.0;
    return;
  }

  model *musub = cov->kappasub[SHAPE_FCT_MEAN];
  double *mu = P(SHAPE_FCT_MEAN);
  if (musub != NULL) {
    FCTN(x, info, musub, v);
  } else for (int i=0; i<vdim; i++) v[i] = ISNAN(mu[i]) ? 0.0 : mu[i];
  // 1.0 notwendig fuer likelihood berechnung;
}

void kappashapefct(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc){
  *nr = i == SHAPE_FCT_MEAN ? SIZE_NOT_DETERMINED : -1; 
  *nc = 1;
}


bool setshapefct(model *cov) {
  model *musub = cov->kappasub[SHAPE_FCT_MEAN]; 
  isotropy_type iso = CONDPREVISO(0); 
  if (!isFixed(iso)) return false;
  set_type(OWN, 0, ShapeType);
  if (musub == NULL) { // derzeit alles coordinatesystems falls trend,
    // da FCTN dies voraussetzt
    set_iso(OWN, 0, PREVISO(0));
    set_xdim(OWN, 0, PREVXDIM(0));
  } else {
    set_iso(OWN, 0, isCartesian(iso) ? CARTESIAN_COORD
	    : isEarth(iso) ? EARTH_COORD
	    : isSpherical(iso) ? SPHERICAL_COORD
	    : ISO_MISMATCH);
    set_xdim(OWN, 0, PREVXDIM(0));
  }
  return true;
}


bool allowedIshapefct(model *cov) {
  model *musub = cov->kappasub[SHAPE_FCT_MEAN]; 
  if (musub == NULL) return allowedItrue(cov);

  bool *I = cov->allowedI;
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  I[CARTESIAN_COORD] = I[EARTH_COORD] = I[SPHERICAL_COORD] = true;
  return false;
}


Types Typeshapefct(Types required, model *cov, isotropy_type required_iso){
  if (!isTrend(required)) return BadType; // 16.1.21 neu
  model *musub = cov->kappasub[SHAPE_FCT_MEAN]; 
  if (musub == NULL) return required;
  return TypeConsistency(required, musub, required_iso);
}


int checkshapefct(model *cov){

  model 
    *calling = cov->calling,
    *musub = cov->kappasub[SHAPE_FCT_MEAN];
  int i,  err,
    vdim = 0, 
    logdim = OWNLOGDIM(0);
    
  if ((musub != NULL) xor PisNULL(SHAPE_FCT_MEAN))
    ERR("exactly one trend argument must be given");

  if (!hasTrendFrame(cov) && !hasAnyEvaluationFrame(cov)) {
     ILLEGAL_FRAME;
  }

  //  PMI(cov);
  if ((cov->matrix_indep_of_x = musub == NULL)) {
    vdim = cov->nrow[SHAPE_FCT_MEAN];
  } else {
    ASSERT_ONESYSTEM;
    assert(musub != NULL);
    if ((err = CHECK(musub, logdim, OWNXDIM(0), ShapeType, XONLY,
		     PREVISO(0), SUBMODEL_DEP, TrendType)) != NOERROR) {
      //    APMI(cov);
      RETURN_ERR(err);
    }
    // PMI(cov);
    // printf("err = %d\n", err);
    vdim = musub->vdim[0];
  }

  VDIM0 = vdim;
  VDIM1 = 1;
  set_iso(OWN, 0, cov->matrix_indep_of_x ? IsotropicOf(OWNISO(0)) 
	  : CoordinateSystemOf(OWNISO(0)));

  RETURN_NOERROR;

}



void rangeshapefct(model  VARIABLE_IS_NOT_USED *cov, range_type *range){
  //P(SHAPE_FCT_MEAN]: mu / mean
  
  range->min[SHAPE_FCT_MEAN] = RF_NEGINF;
  range->max[SHAPE_FCT_MEAN] = RF_INF;
  range->pmin[SHAPE_FCT_MEAN] = -1e10;
  range->pmax[SHAPE_FCT_MEAN] = 1e10;
  range->openmin[SHAPE_FCT_MEAN] = true;
  range->openmax[SHAPE_FCT_MEAN] = true;
}

	

int checkShapefctproc(model *cov) { // auch fuer Shapefctproc
  //  PMI0(cov);
   int err;
  ASSERT_ONESYSTEM;
  // printf("proc %d %d %d %d\n", musub != NULL, PisNULL(SHAPE_FCT_MEAN), cov->sub[0] == NULL, cov->key==NULL);
  if (cov->nsub == 0) {
    if (cov->kappasub[SHAPE_FCT_MEAN] == NULL) ERR("a function must be given");
    assert(PisNULL(SHAPE_FCT_MEAN));
    cov->sub[SHAPE_FCT_MEAN] = cov->kappasub[SHAPE_FCT_MEAN];
    cov->kappasub[SHAPE_FCT_MEAN] = NULL;
    cov->nsub = 1;
  }
  assert(cov->kappasub[SHAPE_FCT_MEAN] == NULL);
  model *musub = cov->sub[SHAPE_FCT_MEAN];
  if ((musub != NULL) xor PisNULL(SHAPE_FCT_MEAN))
    SERR("either a mean 'mu' or an RMmodel must be given");
  if (musub != NULL) {
    int newdim = PREVLOGDIM(0);

    //    PMI(cov);   printf("iso = %s\n", ISO_NAMES[ PREVISO(0)]);
    
    if ((err = CHECK(musub, newdim, newdim, ShapeType, XONLY,
		     PREVISO(0), SUBMODEL_DEP, TrendType)) !=
	NOERROR) RETURN_ERR(err);
    
    setbackward(cov, musub);
    VDIM0 = musub->vdim[0]; 
    VDIM1 = musub->vdim[1];
  } else {
    VDIM0 = cov->nrow[SHAPE_FCT_MEAN];
    VDIM1 = cov->ncol[SHAPE_FCT_MEAN];
  }
  if (VDIM1 != 1) ERR("matrix-valued functions not allowed as trend")

  RETURN_NOERROR;
}
 

int init_Shapefctproc(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){// auch fuer Shapefctproc
  // FRAME_ASSERT_GAUSS;
   model *musub = cov->sub[SHAPE_FCT_MEAN];
   if (VDIM0 != 1) NotProgrammedYet("");
  int err;
  if (musub != NULL && (err = check_fctn_intern(cov)) != NOERROR)
    goto ErrorHandling;
  if ((err = ReturnOwnField(cov)) != NOERROR) goto ErrorHandling;
  if (PL>= PL_STRUCTURE) { PRINTF("\n'%s' is now initialized.\n", NAME(cov));}

 ErrorHandling:  
  cov->simu.active = err == NOERROR;
  RETURN_ERR(err);
}

void do_Shapefctproc(model *cov, gen_storage  VARIABLE_IS_NOT_USED *s){
   model *musub = cov->sub[SHAPE_FCT_MEAN];
  double
    *res = cov->rf;
  assert(res != NULL);
  errorloc_type errorloc_save;
  char *error_location = cov->base->error_location;  
  STRCPY(errorloc_save, error_location);
  SPRINTF(error_location, "%.50s%.50s", errorloc_save, "add trend model");
  if (musub != NULL) Fctn(cov, true, res);
  else {
    int
      vdim = VDIM0,
      vdimtot = Loctotalpoints(cov) * vdim;
    double *mu = (double *) MALLOC(vdim * sizeof(double));
    assert(cov->nrow[SHAPE_FCT_MEAN] == vdim);
    MEMCOPY(mu, P(SHAPE_FCT_MEAN), sizeof(double) * vdim);
    for (int i=0; i<vdimtot; i++) res[i] = mu[i % vdim];
    FREE(mu);
  }
  STRCPY(error_location, errorloc_save);
  return; 
}



bool setshapefctproc(model *cov) {
  model *musub = cov->sub[SHAPE_FCT_MEAN]; 
  isotropy_type iso = CONDPREVISO(0); 
  if (!isFixed(iso)) return false;
  set_type(OWN, 0, ProcessType);
  if (musub == NULL) { // derzeit alles coordinatesystems falls trend,
    // da FCTN dies voraussetzt
    set_iso(OWN, 0, PREVISO(0));
    set_xdim(OWN, 0, PREVXDIM(0));
  } else {
    set_iso(OWN, 0, isCartesian(iso) ? CARTESIAN_COORD
	    : isEarth(iso) ? EARTH_COORD
	    : isSpherical(iso) ? SPHERICAL_COORD
	    : ISO_MISMATCH);
    set_xdim(OWN, 0, PREVXDIM(0));
  }
  return true;
}


bool allowedIshapefctproc(model *cov) {
  model *musub = cov->sub[SHAPE_FCT_MEAN]; 
  if (LocDist(cov) && musub != NULL) return allowedIfalse(cov);

  bool *I = cov->allowedI;
  for (int i=(int) FIRST_ISOUSER; i<=(int) LAST_ISOUSER; I[i++] = false);
  I[CARTESIAN_COORD] = I[EARTH_COORD] = I[SPHERICAL_COORD] = true;
  if (musub == NULL)
    I[ISOTROPIC] = I[EARTH_ISOTROPIC] = I[SPHERICAL_ISOTROPIC] = true;
 return false;
}


Types Typeshapefctproc(Types required, model VARIABLE_IS_NOT_USED*cov,
		    isotropy_type VARIABLE_IS_NOT_USED required_iso){
  return equalsProcess(required) ? required : BadType;
}

