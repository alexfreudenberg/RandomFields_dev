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

# ifndef PrimitivesBi_H
# define PrimitivesBi_H 1

# define UNIT_EPSILON 1E-13
# include "RF.h"
// gauss covariance models
void rangegaussgauss(model *cov, range_type *range);
void gaussgauss(double *x, int *, model *cov, double *v);

void rangegaussGammalike(model *cov, range_type *range);
void gaussGammalike(double *x, int *, model *cov, double *v);

// cauchy covariance models
void rangeCauchyUnif1(model *cov, range_type *range);
void cauchyUnif1(double *x, int *, model *cov, double *v);

void rangeCauchyUnif2(model *cov, range_type *range);
void cauchyUnif2(double *x, int *, model *cov, double *v);

void rangeCauchyUnif3(model *cov, range_type *range);
void cauchyUnif3(double *x, int *, model *cov, double *v);

void rangelatentCauchy1(model *cov, range_type *range);
void latentCauchy1(double *x, int *, model *cov, double *v);

void rangelatentCauchy2(model *cov, range_type *range);
void latentCauchy2(double *x, int *, model *cov, double *v);

void rangelatentCauchy3(model *cov, range_type *range);
void latentCauchy3(double *x, int *, model *cov, double *v);

void rangelatentCauchy4(model *cov, range_type *range);
void latentCauchy4(double *x, int *, model *cov, double *v);

void kappaStrokorbOesting(int i, model *cov, int *nr, int *nc);
int checkStrokorbOesting(model *cov);
void rangeStrokorbOesting(model *cov, range_type *range);
void strokorbOesting(double *x, double *y, int *, model *cov, double *v);

void rangeBiCoshShift(model *cov, range_type *range);
void biCoshShift(double *x, int *, model *cov, double *v);

#endif
