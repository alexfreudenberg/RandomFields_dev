/*
   Authors 
   Felix Reinbott, freinbot@mail.uni-mannheim.de

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

#include <Rmath.h>
#include <R_ext/Lapack.h>

#include "questions.h"
#include "primitive.bivariate.h"
#include "operator.h"
#include "AutoRandomFields.h"
#include "shape.h"
#include "primitive.h"
#include "startGetNset.h"


/* Models involving Gauss */

#define GAUSSGAUSS_NU 0
void gaussgauss(double *x, INFO, model *cov, double *v) {
  // bivariate gauss model
  double x_sum, x_norm_sq,
    nu = P0(GAUSSGAUSS_NU);
  int dim = OWNLOGDIM(0);

  x_norm_sq = 0;
  x_sum = 0;

  // using a for loop instead of the scalar function because it 
  // allows for the computation of both dotproducts.
  for (int i = 0; i < dim; i++) {
    x_sum += x[i];
    x_norm_sq += x[i] * x[i];
  }

  v[0] = EXP(- x_norm_sq);
  v[1] = SQRT(nu / (dim + nu)) * EXP(- x_norm_sq + x_sum * x_sum /(dim + nu));
  v[2] = v[1];
  v[3] = v[0];
}

void rangegaussgauss(model VARIABLE_IS_NOT_USED *cov, range_type *range){
  // nu: in R+
  range->min[GAUSSGAUSS_NU] = 0.0;
  range->max[GAUSSGAUSS_NU] = RF_INF;
  range->pmin[GAUSSGAUSS_NU] = 1e-10;
  range->pmax[GAUSSGAUSS_NU] = 10.0;
  range->openmin[GAUSSGAUSS_NU] = true;
  range->openmax[GAUSSGAUSS_NU] = true;
}

#define GAUSSGAMMALIKE_M 0
void gaussGammalike(double *x, INFO, model *cov, double *v) {
  int dim = OWNLOGDIM(0);
  int m = P0INT(GAUSSGAMMALIKE_M);

  double gamma_factor, inner_factor, 
          sum_value = 0,
         polynomial_factor,
         factorial = 1, 
         x_norm_sq = 0,
         x_sum = 0;

  // calculate the scalar product of x with a ones vector, as well as the squared norm
  for(int i = 0; i < dim; i++) {
    x_norm_sq += x[i] * x[i];
    x_sum += x[i];
  }

  // set the helping factors for the diagonal, power in and outside the sum
  double diag_val = EXP(- x_norm_sq),
         summing_factor = (dim + 1) / 4;


  // calculate the products and factorials for the multiplicative part 
  // before the sum
  for(int i=1; i <= (2 * m); i++) {
    factorial *= i;
  }

  // calculates the sum
  for(int k = 0; k <= m; k++) {
    double x_sum_power,
           fct1 = 1,    // factor for the term with the k-th power
           fct2 = 1,    // for k factorial  
           fct3 = 1;    // for the other factorial term 

    // calculating the power factor
    fct1 = POW(summing_factor, k);

    for(int i = 1; i <= k; i++) {
      fct2 *= i;
    }

    for(int i = 1; i <= (2 * (m - k)); i++) fct3 *= i;

    // caluclate the sum power:
    x_sum_power = POW(x_sum, (2 * (m -  k)));

    // add the kth power term multiplied with the factorial term to the total sum
    sum_value += (fct1 * x_sum_power / (fct2 * fct3));
  }

  polynomial_factor = POW((2 / (dim + 1)), (2 * m));
  gamma_factor = diag_val / gammafn(m + 0.5);
  inner_factor = SQRT(PI / (dim+1)) * factorial * EXP((x_sum * x_sum) / (dim + 1)) * polynomial_factor; 

  // set the diagonal
  v[0] = diag_val;
  v[3] = v[0];

  // multiply all calculated factors together and set off diagonal
  v[1] = gamma_factor * inner_factor * sum_value;
  v[2] = v[1];
}

void rangegaussGammalike(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  // m in |N
  range->min[GAUSSGAMMALIKE_M] = 1;
  // TODO ask for a useful setting to prevent numerical 0's.
  range->max[GAUSSGAMMALIKE_M] = MAXINT;
  range->pmin[GAUSSGAMMALIKE_M] = range->min[GAUSSGAMMALIKE_M];
  range->pmax[GAUSSGAMMALIKE_M] =  10;
  // TODO ask Martin if this is a good idea
  range->openmin[GAUSSGAMMALIKE_M] = false;
  range->openmax[GAUSSGAMMALIKE_M] = false;
}


/* shift models involving cauchy models */
#define CAUCHYUNIF1_EPS 0
#define CAUCHYUNIF1_B 1
void cauchyUnif1(double *x, INFO, model *cov, double *v) {
  double x_sum, x_norm_sq, delta,
         eps = P0(CAUCHYUNIF1_EPS),
         b = P0(CAUCHYUNIF1_B);
  int dim = OWNLOGDIM(0);

  x_norm_sq = 0;
  x_sum = 0;

  // using a for loop instead of the scalar function because it 
  // allows for the computation of both dotproducts.
  for (int i = 0; i < dim; i++) {
    x_sum += x[i];
    x_norm_sq += x[i] * x[i];
  }

  v[0] = 1 / (eps +  x_norm_sq);
  v[3] = v[0];

  delta = SQRT(dim * eps + dim * x_norm_sq - x_sum * x_sum); 
  v[1] = (ATAN((x_sum + dim * b) / delta) - ATAN((x_sum - dim * b) / delta)) / (2 * b * delta);
  v[2] = v[1];
}

void rangeCauchyUnif1(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[CAUCHYUNIF1_EPS] = 0.0;
  range->max[CAUCHYUNIF1_EPS] = RF_INF;
  range->pmin[CAUCHYUNIF1_EPS] = 1e-10;
  range->pmax[CAUCHYUNIF1_EPS] = 10.0;
  range->openmin[CAUCHYUNIF1_EPS] = true;
  range->openmax[CAUCHYUNIF1_EPS] = true;

  range->min[CAUCHYUNIF1_B] = 0.0;
  range->max[CAUCHYUNIF1_B] = RF_INF;
  range->pmin[CAUCHYUNIF1_B] = 1e-10;
  range->pmax[CAUCHYUNIF1_B] = 10.0;
  range->openmin[CAUCHYUNIF1_B] = true;
  range->openmax[CAUCHYUNIF1_B] = true;
}


#define CAUCHYUNIF2_EPS 0
#define CAUCHYUNIF2_B 1
void cauchyUnif2(double *x, INFO, model *cov, double *v) {
  double x_sum, x_norm_sq, root_val, fract_const, val_1, val_2,
         eps = P0(CAUCHYUNIF2_EPS),
         b = P0(CAUCHYUNIF2_B);
  int dim = OWNLOGDIM(0);

  x_norm_sq = 0;
  x_sum = 0;

  // using a for loop instead of the scalar function because it 
  // allows for the computation of both dotproducts.
  for (int i = 0; i < dim; i++) {
    x_sum += x[i];
    x_norm_sq += x[i] * x[i];
  }

  root_val = SQRT(eps +  x_norm_sq);
  v[0] = 1/ root_val;
  v[3] = v[0];

  // off diagonal
  fract_const = SQRT(dim * (eps + x_norm_sq) - x_sum * x_sum);
  val_1 = (dim * b + x_sum) / fract_const;
  val_2 = (-dim * b + x_sum) / fract_const;

  // using the naive formula for arsinh 
  v[1] = (LOG(val_1 + SQRT(val_1 * val_1 + 1)) - LOG(val_2 + SQRT(val_2 * val_2 + 1))) / (2 * b * SQRT(dim));
  v[2] = v[1];
}

void rangeCauchyUnif2(model VARIABLE_IS_NOT_USED  *cov, range_type *range) {
  range->min[CAUCHYUNIF2_EPS] = 0.0;
  range->max[CAUCHYUNIF2_EPS] = RF_INF;
  range->pmin[CAUCHYUNIF2_EPS] = 1e-10;
  range->pmax[CAUCHYUNIF2_EPS] = 10.0;
  range->openmin[CAUCHYUNIF2_EPS] = true;
  range->openmax[CAUCHYUNIF2_EPS] = true;

  range->min[CAUCHYUNIF2_B] = 0.0;
  range->max[CAUCHYUNIF2_B] = RF_INF;
  range->pmin[CAUCHYUNIF2_B] = 1e-10;
  range->pmax[CAUCHYUNIF2_B] = 10.0;
  range->openmin[CAUCHYUNIF2_B] = true;
  range->openmax[CAUCHYUNIF2_B] = true;
}

#define CAUCHYUNIF3_EPS 0
#define CAUCHYUNIF3_B 1
void cauchyUnif3(double *x, INFO, model *cov, double *v) {
  double x_sum, x_norm_sq, root_val,
         eps = P0(CAUCHYUNIF3_EPS),
         b = P0(CAUCHYUNIF3_B);
  int dim = OWNLOGDIM(0);

  x_norm_sq = 0;
  x_sum = 0;

  // using a for loop instead of the scalar function because it 
  // allows for the computation of both dotproducts.
  for (int i = 0; i < dim; i++) {
    x_sum += x[i];
    x_norm_sq += x[i] * x[i];
  }

  root_val = SQRT(eps +  x_norm_sq);
  v[0] = 1/ (root_val * root_val * root_val);
  v[3] = v[0];

  // talk with martin about that
  double c1, c2, c3, counter1, counter2;

  counter1 = (dim * b + x_sum);
  counter2 = (- dim * b + x_sum);
  c1 = 1 / (dim * (eps + x_norm_sq) - x_sum * x_sum);
  c2 = 1 / SQRT(eps + x_norm_sq + 2 * x_sum * b + dim * b * b);
  c3 = 1 / SQRT(eps + x_norm_sq - 2 * x_sum * b + dim * b * b);
  v[1] = 1 / (2 * b) * ((counter1 * c1) * c2 - (counter2 * c1) * c3);
  v[2] = v[1];
}


void rangeCauchyUnif3(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[CAUCHYUNIF3_EPS] = 0.0;
  range->max[CAUCHYUNIF3_EPS] = RF_INF;
  range->pmin[CAUCHYUNIF3_EPS] = 1e-10;
  range->pmax[CAUCHYUNIF3_EPS] = 10.0;
  range->openmin[CAUCHYUNIF3_EPS] = true;
  range->openmax[CAUCHYUNIF3_EPS] = true;

  range->min[CAUCHYUNIF3_B] = 0.0;
  range->max[CAUCHYUNIF3_B] = RF_INF;
  range->pmin[CAUCHYUNIF3_B] = 1e-10;
  range->pmax[CAUCHYUNIF3_B] = 10.0;
  range->openmin[CAUCHYUNIF3_B] = true;
  range->openmax[CAUCHYUNIF3_B] = true;
}

# define LATENTCAUCHY1_A 0
# define LATENTCAUCHY1_B 1
void latentCauchy1(double *x, INFO, model *cov, double *v) {
  // declare variables
  int dim = OWNLOGDIM(0);
  double x_norm_sq = 0;
  double a = P0(LATENTCAUCHY1_A);
  double b = P0(LATENTCAUCHY1_B);

  for(int i = 0; i < dim; i++) {
    x_norm_sq += x[i] * x[i];
  }

  v[0] = 1 / (1 + x_norm_sq);
  v[3] = v[0];

  v[1] = 1 / (b - a) * SQRT(v[0]) * (ATAN(b * (SQRT(v[0]))) - ATAN(a * SQRT(v[0])));
  v[2] = v[1];
}

void rangelatentCauchy1(model VARIABLE_IS_NOT_USED  *cov, range_type *range) {
  range->min[LATENTCAUCHY1_A] = -RF_INF;
  range->max[LATENTCAUCHY1_A] = RF_INF;
  range->pmin[LATENTCAUCHY1_A] = -1000;
  range->pmax[LATENTCAUCHY1_A] = 1000;
  range->openmin[LATENTCAUCHY1_A] = true;
  range->openmax[LATENTCAUCHY1_A] = true;

  range->min[LATENTCAUCHY1_B] = -RF_INF;
  range->max[LATENTCAUCHY1_B] = RF_INF;
  range->pmin[LATENTCAUCHY1_B] = -1000;
  range->pmax[LATENTCAUCHY1_B] = 1000;
  range->openmin[LATENTCAUCHY1_B] = true;
  range->openmax[LATENTCAUCHY1_B] = true;
}

# define LATENTCAUCHY2_A 0
# define LATENTCAUCHY2_B 1
void latentCauchy2(double *x, INFO, model *cov, double *v) {
  // declare variables
  int dim = OWNLOGDIM(0);
  double x_norm_sq = 0, 
         a = P0(LATENTCAUCHY2_A), 
         b = P0(LATENTCAUCHY2_B);


  for(int i = 0; i < dim; i++) {
    x_norm_sq += x[i] * x[i];
  }

  v[0] = 1 / SQRT(1 + x_norm_sq);
  v[3] = v[0];

  v[1] = 1 / (b - a) * LOG((b + SQRT(1 + b*b + x_norm_sq)) / ( a + SQRT(1 + a*a + x_norm_sq)));
  v[2] = v[1];
}

void rangelatentCauchy2(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[LATENTCAUCHY2_A] = -RF_INF;
  range->max[LATENTCAUCHY2_A] = RF_INF;
  range->pmin[LATENTCAUCHY2_A] = -1000;
  range->pmax[LATENTCAUCHY2_A] = 1000;
  range->openmin[LATENTCAUCHY2_A] = true;
  range->openmax[LATENTCAUCHY2_A] = true;

  range->min[LATENTCAUCHY2_B] = -RF_INF;
  range->max[LATENTCAUCHY2_B] = RF_INF;
  range->pmin[LATENTCAUCHY2_B] = -1000;
  range->pmax[LATENTCAUCHY2_B] = 1000;
  range->openmin[LATENTCAUCHY2_B] = true;
  range->openmax[LATENTCAUCHY2_B] = true;
}

# define LATENTCAUCHY3_A 0
# define LATENTCAUCHY3_B 1
void latentCauchy3(double *x, INFO, model *cov, double *v) {
  // declare variables
  int dim = OWNLOGDIM(0);
  double x_norm_sq = 0, 
         intermed,
         a = P0(LATENTCAUCHY3_A), 
         b = P0(LATENTCAUCHY3_B);


  for(int i = 0; i < dim; i++) {
    x_norm_sq += x[i] * x[i];
  }

  intermed = 1 / SQRT(1 + x_norm_sq);
  v[0] = intermed * intermed * intermed;
  v[3] = v[0];

  v[1] = 1 / ((b - a) * (1 + x_norm_sq)) *
    (b / (SQRT(1 + b*b + x_norm_sq)) - a / (SQRT(1 + a*a + x_norm_sq))); 
  v[2] = v[1];
}

void rangelatentCauchy3(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  // 
  range->min[LATENTCAUCHY3_A] = -RF_INF;
  range->max[LATENTCAUCHY3_A] = RF_INF;
  range->pmin[LATENTCAUCHY3_A] = -1000;
  range->pmax[LATENTCAUCHY3_A] = 1000;
  range->openmin[LATENTCAUCHY3_A] = true;
  range->openmax[LATENTCAUCHY3_A] = true;

  range->min[LATENTCAUCHY3_B] = -RF_INF;
  range->max[LATENTCAUCHY3_B] = RF_INF;
  range->pmin[LATENTCAUCHY3_B] = -1000;
  range->pmax[LATENTCAUCHY3_B] = 1000;
  range->openmin[LATENTCAUCHY3_B] = true;
  range->openmax[LATENTCAUCHY3_B] = true;
}

# define LATENTCAUCHY4_A 0
# define LATENTCAUCHY4_B 1
# define LATENTCAUCHY4_GAMMA 2
void latentCauchy4(double *x, INFO, model *cov, double *v) {
  // declare variables
  int dim = OWNLOGDIM(0);
  double x_norm_sq = 0, 
         gamma = P0(LATENTCAUCHY4_GAMMA),
         a = P0(LATENTCAUCHY4_A), 
         b = P0(LATENTCAUCHY4_B);


  for(int i = 0; i < dim; i++) {
    x_norm_sq += x[i] * x[i];
  }

  v[0] = 1 / POW(1 + x_norm_sq, gamma);
  v[3] = v[0];

  v[1] = 1 / (b*b - a*a) * (gamma - 1) * (1 / POW(1 + a*a + x_norm_sq, gamma - 1) - 1 / POW(1 + b*b + x_norm_sq, gamma - 1));
  v[2] = v[1];
}

void rangelatentCauchy4(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  // 
  range->min[LATENTCAUCHY4_A] = 0;
  range->max[LATENTCAUCHY4_A] = RF_INF;
  range->pmin[LATENTCAUCHY4_A] = 0;
  range->pmax[LATENTCAUCHY4_A] = 1000;
  range->openmin[LATENTCAUCHY4_A] = true;
  range->openmax[LATENTCAUCHY4_A] = true;

  range->min[LATENTCAUCHY4_B] = 0;
  range->max[LATENTCAUCHY4_B] = RF_INF;
  range->pmin[LATENTCAUCHY4_B] = 0;
  range->pmax[LATENTCAUCHY4_B] = 1000;
  range->openmin[LATENTCAUCHY4_B] = true;
  range->openmax[LATENTCAUCHY4_B] = true;
  
  range->min[LATENTCAUCHY4_GAMMA] = 0;
  range->max[LATENTCAUCHY4_GAMMA] = RF_INF;
  range->pmin[LATENTCAUCHY4_GAMMA] = 0;
  range->pmax[LATENTCAUCHY4_GAMMA] = 1000;
  range->openmin[LATENTCAUCHY4_GAMMA] = true;
  range->openmax[LATENTCAUCHY4_GAMMA] = true;
}

/* Miscellaneous models */ 
# define STROES_X0 0
# define STROES_Y0 1
# define STROES_PHI 2
# define STROES_PSI 3
void strokorbOesting(double *x, int *info, double *y, model *cov, double *v) {
  // TODO make the TBM work... (sorry martin...)

  int dim = OWNLOGDIM(0);

  // starting off with setting X,Y as exponential models
  double *X0 = P(STROES_X0), 
         *Y0 = P(STROES_Y0),
         phi_cov, psi_cov, 
         dx0y0 = 0, 
         dxy0 = 0, 
         dyx0 = 0, 
         dxx0 = 0, 
         dyy0 = 0, 
         dxy = 0;

  model *phi = cov->sub[0], 
        *psi = cov->sub[1];

  for(int i = 0; i < dim; i++) {
    dx0y0 += (X0[i] - Y0[i]) * (X0[i] - Y0[i]);
    dxy0 += (x[i] - Y0[i]) * (x[i] - Y0[i]);
    dyx0 += (y[i] - X0[i]) * (y[i] - X0[i]);
    dyy0 += (y[i] - Y0[i]) * (y[i] - Y0[i]);
    dxx0 += (x[i] - X0[i]) * (x[i] - X0[i]);
    dxy += (x[i] - y[i]) * (x[i] - y[i]);
  }


  dx0y0 = SQRT(dx0y0);
  dxy0 = SQRT(dxy0);
  dyx0 = SQRT(dyx0);
  dxx0 = SQRT(dxx0);
  dyy0 = SQRT(dyy0);
  dxy = SQRT(dxy);

  COV(&dxy, info, phi, &phi_cov);
  COV(&dx0y0, info, psi, &psi_cov);
  v[0] = phi_cov * psi_cov;

  COV(&dxy0, info, phi, &phi_cov);
  COV(&dyx0, info, psi, &psi_cov);
  v[1] = phi_cov * psi_cov;

  COV(&dyy0, info, phi, &phi_cov);
  COV(&dxx0, info, psi, &psi_cov);
  v[2] = phi_cov * psi_cov;

  COV(&dx0y0, info, phi, &phi_cov);
  COV(&dxy, info, psi, &psi_cov);
  v[3] = phi_cov * psi_cov;

}

void kappaStrokorbOesting(int i, model *cov, int *nr, int *nc) {
  *nc = 1;
  *nr = i < DefList[COVNR].kappas ? OWNLOGDIM(0) : 1;
} 

int checkStrokorbOesting(model *cov) {
  int err, dim = OWNLOGDIM(0);
  model *phi = cov ->sub[0], *psi = cov->sub[1];

  if ((err = CHECK(phi, dim, 1, PosDefType, XONLY, ISOTROPIC, 
          1, EvaluationType)) != NOERROR) {
    RETURN_ERR(err); 
  }

  if ((err = CHECK(psi, dim, 1, PosDefType, XONLY, ISOTROPIC,
          1, EvaluationType)) != NOERROR) {
    RETURN_ERR(err); 
  }

  RETURN_NOERROR;
}

void rangeStrokorbOesting(model VARIABLE_IS_NOT_USED *cov, range_type *range) {
  range->min[STROES_X0] = -1000;
  range->max[STROES_X0] = RF_INF;
  range->pmin[STROES_X0] = -1000;
  range->pmax[STROES_X0] = 1000;
  range->openmin[STROES_X0] = true;
  range->openmax[STROES_X0] = true;

  range->min[STROES_Y0] = -1000;
  range->max[STROES_Y0] = RF_INF;
  range->pmin[STROES_Y0] = -1000;
  range->pmax[STROES_Y0] = 1000;
  range->openmin[STROES_Y0] = true;
  range->openmax[STROES_Y0] = true;

}


void includeAsymmetricModels() {
  pref_type pbivariate = {5, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 5};
  //                      CE CO CI TBM Sp di sq Ma av n mpp Hy spf any
  IncludePrim("gaussgauss", PosDefType, 1, XONLY, SYMMETRIC, checkOK, rangegaussgauss, pbivariate, 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
  kappanames("nu", REALSXP);
  addCov(gaussgauss);
 
  IncludePrim("gaussGammalike", PosDefType, 1, XONLY, SYMMETRIC, checkOK,
	      rangegaussGammalike, pbivariate ,
	      2/*this is hopefully the parameter that controls the dim*/,
	      INFDIM, (ext_bool) false, NOT_MONOTONE);
  kappanames("m", INTSXP); // TODO ask or find version for int
  addCov(gaussGammalike); 

  
    IncludePrim("cauchyUnif1", PosDefType, 2, XONLY, SYMMETRIC, checkOK, rangeCauchyUnif1, pbivariate , 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
    kappanames("eps", REALSXP, "b", REALSXP);
    addCov(cauchyUnif1);


     IncludePrim("cauchyUnif2", PosDefType, 2, XONLY, SYMMETRIC, checkOK, rangeCauchyUnif2, pbivariate , 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
    kappanames("eps", REALSXP, "b", REALSXP);
    addCov(cauchyUnif2);


     IncludePrim("cauchyUnif3", PosDefType, 2, XONLY, SYMMETRIC, checkOK, rangeCauchyUnif3, pbivariate , 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
    kappanames("eps", REALSXP, "b", REALSXP);
    addCov(cauchyUnif3);

     IncludePrim("latentCauchy1", PosDefType, 2, XONLY, SYMMETRIC, checkOK, rangelatentCauchy1, pbivariate ,2, INFDIM, (ext_bool) false, NOT_MONOTONE);
    kappanames("a", REALSXP, "b", REALSXP);
    addCov(latentCauchy1);

     IncludePrim("latentCauchy2", PosDefType, 2, XONLY, SYMMETRIC, checkOK, rangelatentCauchy2, pbivariate , 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
    kappanames("a", REALSXP, "b", REALSXP);
    addCov(latentCauchy2);

    IncludePrim("latentCauchy3", PosDefType, 2, XONLY, SYMMETRIC, checkOK, rangelatentCauchy3, pbivariate , 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
    kappanames("a", REALSXP, "b", REALSXP);
    addCov(latentCauchy3);

   IncludePrim("latentCauchy4", PosDefType, 3, XONLY, SYMMETRIC, checkOK, rangelatentCauchy4, pbivariate , 2, INFDIM, (ext_bool) false, NOT_MONOTONE);
    kappanames("a", REALSXP, "b", REALSXP, "gamma", REALSXP);
    addCov(latentCauchy4);

}
