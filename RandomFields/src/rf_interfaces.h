
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

#ifndef RF_INTERFACE_H
#define RF_INTERFACE_H 1


model *InitIntern(int cR, SEXP Model, SEXP x, bool NA_OK,
		  raw_type rawConcerns);

void simulate(double *x, int*, model *cov, double *v);
int check_simulate(model *cov); 
void range_simulate(model *cov, range_type *range);  
int struct_simulate(model *cov, model **newmodel);
int init_simulate(model *cov, gen_storage *S);
// void do_simulate(model *cov, gen_storage *S);
void range_simulate(model VARIABLE_IS_NOT_USED *cov, range_type* range);

void density(double *x, int*, model *cov, double *v);
int check_density(model *cov); 
void range_density(model *cov, range_type *range);  
int struct_density(model *cov, model **newmodel);
int init_density(model *cov, gen_storage *S);
// void do_density(model *cov, gen_storage *S);
void range_density(model VARIABLE_IS_NOT_USED *cov, range_type* range);



#define LIKELIHOOD_DATA 0
#define LIKELIHOOD_NA_VAR 1
#define LIKELIHOOD_BETASSEPARATE 2
#define LIKELIHOOD_IGNORETREND 3
#define LIKELIHOOD_STANDARDIZED_L 4
#define LIKELIHOOD_LAST LIKELIHOOD_STANDARDIZED_L
#define LIKELI_NA_INTEGER false
#define LIKELI_EXCLUDE_TREND true
void kappalikelihood(int i, model VARIABLE_IS_NOT_USED *cov, 
		     int *nr, int *nc);
void likelihood(double *data, int*, model *cov, double *v);
int check_likelihood(model *cov);
int struct_likelihood(model *cov, model **newmodel);
void range_likelihood(model *cov, range_type* range);

void linearpart(double *data, int *, model *cov, double *v);
int check_linearpart(model *cov);
int struct_linearpart(model *cov, model **newmodel);
//void range_linearpart(model *cov, range_type* range);

#define PREDICT_DATA 0
#define PREDICT_NA_VAR 1
#define PREDICT_BETASSEPARATE 2
#define PREDICT_IGNORETREND 3
#define PREDICT_STANDARDIZED_L 4
#define PREDICT_GIVEN 5
//#define PREDICT_PREDICTIDX 6
#define PREDICT_CONDITIONING 0
#define PREDICT_PREDICT 1

void kappapredict(int i, model VARIABLE_IS_NOT_USED *cov, int *nr, int *nc);
void predict(double VARIABLE_IS_NOT_USED *x, int *, model *cov, double *v);
int check_predict(model *cov);
int struct_predict(model *cov, model VARIABLE_IS_NOT_USED  **newmodel);
void range_predict(model VARIABLE_IS_NOT_USED *predict, range_type* range);

void Cov(double *x, int *,  model *cov, double *v);
int check_cov(model *cov);
int struct_cov(model *cov, model **newmodel);
int init_cov(model *cov, gen_storage *s);


void FctnIntern(model *cov, model *covVdim, model *sub, bool ignore_y,
		double *v);
void FctnExtern(model *cov, model *sub, bool ignore_y, double *v);
void Fctn(double *x, int *, model *cov, double *v);
void Fctn(model *cov, bool ignore_y, double *v);
int check_fctn_intern(model *cov);
int check_fctn(model *cov);

void CovMatrix(double *x, int *, model *cov, double *v);
int check_covmatrix(model *cov) ;

void EvalDistr(double *x, int*, model *cov, double *v);
void kappa_EvalDistr(int i, model *cov, int *nr, int *nc);
int check_EvalDistr(model *cov); 
void range_EvalDistr(model *cov, range_type *range);  
int struct_EvalDistr(model *cov, model **newmodel);
int init_EvalDistr(model *cov, gen_storage *S);
// void do_EvalDistr(model *cov, gen_storage *S);

void RFget(double *x, int*, model *cov, double *v);
int SearchParam(model *cov, get_storage *s, model_storage *STOMODEL) ;
int check_RFget(model *cov) ;
void range_RFget(model *cov, range_type* range);
int struct_RFget(model *cov, model **newmodel);


#define PSEUDO_ALPHA 0
void Pseudovariogram(double *x, int*, model *cov, double *v) ;
int check_pseudovario(model *cov);
void range_pseudovario(model *cov, range_type* range);
void Pseudomadogram(double *x, int *, model *cov, double *v) ;
int check_pseudomado(model *cov);
void range_pseudomado(model *cov, range_type* range);

void Variogram(double *x, int*, model *cov, double *v) ;
int check_vario(model *cov);
int struct_variogram(model *cov, model **newmodel);

void Dummy(double *x, int*, model *cov, double *v);
int check_dummy(model *cov);
int struct_dummy(model *cov, model **newmodel);


//-----------------------------------------------------------------
// unsorted

#endif
