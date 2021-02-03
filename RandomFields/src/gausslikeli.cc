
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Definition of correlation functions and derivatives of hypermodels

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

#include <Rmath.h>
#include "questions.h"
#include "primitive.h"
#include "operator.h"
#include "Processes.h"
#include "variogramAndCo.h"
#include "rf_interfaces.h"
#include "shape.h"
#include "QMath.h"




// *sub = cov->key == NULL ? cov->sub[0] : cov->key;
#define SCALAR_P(A,B,C) Ext_scalarX(A,B,C, SCALAR_AVX)

#define START_GAUSS_INTERFACE						\
   int cR = INTEGER(model_reg)[0];					\
  if (cR < 0 || cR > MODEL_MAX) BUG;					\
  set_currentRegister(cR);						\
  model *cov = KEY()[cR];						\
  assert(cov != NULL && equalsnowInterface(cov)) 			\
  cov = cov->key != NULL ? cov->key : cov->sub[0];			\
  if (COVNR != GAUSSPROC)						\
    ERR("register not initialised as Gaussian likelihood");		\
  model *calling = cov->calling;					\
  if (calling == NULL ||						\
      (CALLINGNR != LIKELIHOOD_CALL && CALLINGNR != LINEARPART_CALL && \
       CALLINGNR != PREDICT_CALL)					\
      ) BUG;								\
  getStorage(L, likelihood);						\
  cov->base->set = 0


int matrixcopyNA(double *dest, double *src, double *cond,
		 int rows, int cols, int append_cond) {
  int k = 0;
  double *ptsrc = src;
  for (int j=0; j<cols; j++) {
    for (int i=0; i<rows; i++, ptsrc++) {
      if (!ISNAN(cond[i])) dest[k++] = *ptsrc;
    }
  }
  
  double *C = cond;
  int m = 0;
  for (int j=0; j<append_cond; j++) {
    for (int i=0; i<rows; i++, m++) {
      if (!ISNAN(C[m])) dest[k++] = C[m];
    }
  }

  if (k == 0) RFERROR("one of the data set seems to consist of NAs only");
  int notnas = k == 0 ? 0 : k / (append_cond + cols); 
  
  assert(notnas * (append_cond + cols) == k);
  return notnas;
}


void SqMatrixcopyNA(double *dest, double *src, double *cond, int rows) {
  int k = 0;
  for (int j=0; j<rows; j++) {
    if (!ISNAN(cond[j])) {
      double *ptsrc = src + j * rows;     
      for (int i=0; i<rows; i++, ptsrc++) {
	if (!ISNAN(cond[i]))  dest[k++] = *ptsrc;
      }
    }
  }
}


 
void boxcox_inverse(double boxcox[], int vdim, double *res, int pts, 
		    int repet) {
  assert(boxcox != NULL);
  if (boxcox[0] == RF_INF) return;
  // (y * L + 1)^{1/L} - mu
  int totvdim = vdim * pts;
  for (int v=0; v<vdim; v++) {
    double lambda = boxcox[2 * v];
    double *p = res + v * pts,
      invlambda = 1.0 / lambda,
      mu = boxcox[2 * v + 1];
    if (!ISNA(lambda)) {
      if (lambda == 1.0 && mu == 1.0) continue;
      if (FABS(lambda) < 1e-20) {	
	for (int r=0; r<repet; r++, p+=totvdim) {
	  for (Long i=0; i<pts; i++) {
	    p[i] = EXP(p[i]) - mu; 
	  }
	}
      }
    } else {      
      for (int r=0; r<repet; r++, p+=totvdim) {
	for (int i=0; i<pts; i++) {
	  double dummy = lambda * p[i] + 1.0;
	  if ((dummy < 0 && lambda != CEIL(lambda)) || 
	      (dummy == 0 && invlambda <= 0) ) {
	    RFERROR("value(s) in the inverse Box-Cox transformation not positive");
	  }
	  p[i] = POW(dummy, invlambda) - mu;
	}
      }
    }
  }
}

void boxcox_trafo(double boxcox[], int vdim, double *res, Long pts, int repet){	
  assert(boxcox != NULL);
  if (boxcox[0] == RF_INF) return;
  // box cox: (x + mu)^L - 1) / L
  int totvdim = vdim * pts;
  for (int v=0; v<vdim; v++) {
    double lambda = boxcox[2 * v];
    double *p = res + v * pts,
      invlambda = 1.0 / lambda,
	mu = boxcox[2 * v + 1];
    if (!ISNA(lambda)) {
      if (lambda == 1.0 && mu == 1.0) continue;
      if (FABS(lambda) < 1e-20) {	
	// to do taylor
	for (int r=0; r<repet; r++, p+=totvdim) {
	  for (Long i=0; i<pts; i++) {
	    double dummy = res[i] + mu;
	    if (dummy < 0 || (dummy == 0 && lambda <= 0) ) {
	      RFERROR2("%ld-th value in the Box-Cox transformation not positive (%e)", i, dummy); // PRINTF
	    }
	    res[i] = LOG(dummy);			
	  }
	}
      }
    } else {
      for (int r=0; r<repet; r++, p+=totvdim) {
	for (Long i=0; i<pts; i++) {
	  //printf("%d %f\n", i, res[i]);
	  double dummy = res[i] + mu;
	  if ((dummy < 0 && lambda != CEIL(lambda)) 
	      || (dummy == 0 && lambda <= 0) ) {
	    RFERROR2("%ld-th value in the Box-Cox transformation not positive (%e).", i, dummy);// PRINTF
	  }
	  res[i] = (POW(dummy, lambda) - 1.0) * invlambda; 
	}
      }
    }
  }
}


SEXP BoxCox_trafo(SEXP boxcox, SEXP res, SEXP Vdim, SEXP inverse){	
  int vdim = INTEGER(Vdim)[0],
    repet = isVector(res) ? 1 : ncols(res),
    pts = isVector(res) ? length(res) / vdim : nrows(res);
  if (length(res) != pts * vdim * repet)
    RFERROR("multi-dimensionality incorrect in Box-Cox transformation");
  if (length(boxcox) < 2 * vdim) RFERROR("too few entries in 'boxcox'");
  if (LOGICAL(inverse)[0])
    boxcox_inverse(REAL(boxcox), vdim, REAL(res), pts, repet);
  else 
    boxcox_trafo(REAL(boxcox), vdim, REAL(res), pts, repet);
  return R_NilValue;
}


SEXP set_boxcox(SEXP boxcox, SEXP Reg) {
  globalparam *global = &GLOBAL;
  if (length(Reg) > 0) {
    int reg = INTEGER(Reg)[0];
    set_currentRegister(reg);
    KEY_type *KT = KEYT();
    if (KT == NULL) BUG;
    global = &(KT->global);
  }
  double *bc = REAL(boxcox);
  int len = length(boxcox);
  for (int i=0; i<len; i++) global->gauss.boxcox[i] = bc[i];
  global->gauss.loggauss = false;
  return R_NilValue;
}


SEXP get_boxcox(SEXP Reg) {
  globalparam *global = &GLOBAL;
  if (length(Reg) > 0) {
    int reg = INTEGER(Reg)[0];
    set_currentRegister(reg);
    KEY_type *KT = KEYT();
    if (KT == NULL) BUG;
    global = &(KT->global);
  }
  SEXP ans;
  int len = 2 * MAXBOXCOXVDIM;
  PROTECT(ans =allocVector(REALSXP, len));
  if (global->gauss.loggauss) {
    for (int i=0; i<len; REAL(ans)[i++] = 0.0);
  } else  for (int i=0; i<len; i++) REAL(ans)[i] = global->gauss.boxcox[i];
  UNPROTECT(1);
  return ans;
}


SEXP gauss_linearpart(SEXP model_reg, SEXP Set){
  START_GAUSS_INTERFACE;
#define ll 3
  const char *names[ll] = {"Y", "X", "vdim"};
  SEXP ans, namevec, X, Y;
  int 
    unp = 4, // ans, namevec, X, Y
    sets = L->sets, 
    element  = INTEGER(Set)[0],
    vdim = VDIM0,
    betatot = L->cum_n_betas[L->fixedtrends];

  if (element > 0  && element > sets) ERR("set number out of range");

  PROTECT(ans = allocVector(VECSXP, ll));
  PROTECT(namevec = allocVector(STRSXP, ll));
  for (int k=0; k<ll; k++) SET_STRING_ELT(namevec, k, mkChar(names[k]));
  
  if (element <= 0) {
    PROTECT(Y = allocVector(VECSXP, sets)); 
    PROTECT(X = allocVector(VECSXP, sets)); 
    for (int set = 0; set < sets; set++) {
      cov->base->set =set;
      int totptsvdim = LoctotalpointsY(cov) * vdim;
      SEXP partY , partX;
      PROTECT(partY = allocVector(REALSXP, totptsvdim));
      MEMCOPY(REAL(partY), L->YhatWithoutNA[set], 
	      totptsvdim * sizeof(double));
      SET_VECTOR_ELT(Y, set, partY);
      UNPROTECT(1);
      
      if (L->fixedtrends) {
	PROTECT(partX = allocMatrix(REALSXP, totptsvdim, betatot));
	MEMCOPY(REAL(partX), L->X[set], 
		totptsvdim * betatot * sizeof(double));
	SET_VECTOR_ELT(X, set, partX);
	UNPROTECT(1);
      } else {
	SET_VECTOR_ELT(ans, set, R_NilValue);    
      }
    }
  } else {
    cov->base->set = element - 1;
    int totptsvdim = LoctotalpointsY(cov) * vdim;
   
    PROTECT(Y = allocVector(REALSXP, totptsvdim));
    MEMCOPY(REAL(Y), 
	    L->YhatWithoutNA[cov->base->set], totptsvdim * sizeof(double));
  
    if (L->fixedtrends) {
      PROTECT(X = allocMatrix(REALSXP, totptsvdim, betatot));
      MEMCOPY(REAL(X), L->X[cov->base->set],
	       totptsvdim * betatot * sizeof(double));
    } else {
      X = R_NilValue; 
      unp--;
    }
  }
  
  int k=0;
  SET_VECTOR_ELT(ans, k++, Y);    
  SET_VECTOR_ELT(ans, k++, X);    
  SET_VECTOR_ELT(ans, k++, ScalarInteger(vdim));
  assert(k == ll);

  setAttrib(ans, R_NamesSymbol, namevec);
  UNPROTECT(unp);

  cov->base->set = 0;
  return ans;
}

 

#define LOGLI_RESIDUALS 0
#define LOGLI_WHOLETREND 1
void get_logli_residuals(model *cov, double *work, double *ans,
			 int modus) {
  // do NOT set cov->base->set = 0;
  assert(isnowProcess(cov));
  getStorage(L ,   likelihood);
  assert(L != NULL);
  listoftype *datasets = L->datasets;
  assert(datasets != NULL && ans != NULL);

  int
    vdim = VDIM0,
    betas = L->cum_n_betas[L->fixedtrends],
    ncol = NCOL_OUT_OF(datasets),
    repet = ncol / vdim,
    nrow = NROW_OUT_OF(datasets),
    ndata = nrow * ncol,
    totptsvdim = nrow * vdim;
  double *X = L->X[cov->base->set],
    *pres = ans,
    *data = SET_OUT_OF(datasets);

  if (modus == LOGLI_RESIDUALS) {
    MEMCOPY(pres, data, ndata * sizeof(double));
    boxcox_trafo(P(GAUSS_BOXCOX), vdim, pres, nrow, repet);
  } else {
    // BUG; printf("what's that good for?");
    for (int i=0; i<ndata; pres[i++] = 0.0);
  }

  //   for (int i=0; i<ndata; i++) printf("%10g ", pres[i]);  BUG;
  //  printf("\n%d %d %d vdim=%d\n", L->ignore_trend, L->dettrends, L->fixedtrends, vdim);

  if (L->ignore_trend) {
     return;
  }

  bool delete_work;
  if ((delete_work = work == NULL)) {
    work = (double *) MALLOC(nrow * vdim * sizeof(double));
  }

  double *betavec = L->betavec;
  int r, z;
  
  if (L->dettrends) {
    for (int i=0; i<L->dettrends; i++) {      
      if (L->nas_det[i]) {// non-linear parts without linear part
	FctnIntern(cov, L->det_effect[i], L->det_effect[i], false, work);
	for (z=r=0; r<repet; r++)
	  for (int k=0; k<totptsvdim; k++) pres[z++] -= work[k];
      }
    }
    int set = cov->base->set;
    for (z=r=0; r<repet; r++)
      for (int k=0; k<totptsvdim; k++)
	pres[z++] -= L->YhatWithoutNA[set][k];
  }
  
  if (L->fixedtrends) {
    for (r=0; r<repet; r++, betavec += betas) {
      if (r==0 || L->betas_separate) {
	for (int i=0; i<totptsvdim; work[i++] = 0.0);
	//  for (int j=0; j<totptsvdim;j ++) printf("j=%d p %10g d %10g w %10g\n ", j, pres[j], data[j], work[j]);
	for (int i=0; i<betas; i++) {
	  double beta = betavec[i];
	  //	    printf("i=%d beta = %10g %d\n", i, beta, nrow);
	  for (int j=0; j<nrow; j++, X++) {	
	    //	      printf("i=%d %d %10g %10g\n", i, betas, beta, *X);
	    work[j] += *X * beta;
	  }
	}
      }
      
      for (int j=0; j<nrow; j++, data++, pres++) {
	  *pres -= work[j];
	  //    printf("j=%d p %10g d %10g w %10g\n ", j, *pres, *data, work[j]);
      }
    }
  }
  
  // for (int j=0; j<nrow; j++, data++, pres++) printf("%10g %10g\n", *pres, *data);
  if (modus != LOGLI_RESIDUALS) {
    //BUG; printf("whats that good for?");
    for (int i=0; i<ndata; i++) pres[i] = -pres[i];
  }
  
  if (delete_work) FREE(work);
}


SEXP get_logli_residuals(model *cov, int modus) {
  cov->base->set = 0;
  getStorage(L ,   likelihood);
  listoftype *datasets = L->datasets;
  int 
    vdim = VDIM0,
    sets = LocSets(cov);
  bool matrix = false;
  SEXP all_res, res;

  assert(!L->betas_separate || sets == 1);

  Long max = 0;
  for (int set = 0; set < sets; set++){
    cov->base->set = set;
    Long
      nrow = NROW_OUT_OF(datasets),
      totptsvdim = nrow * vdim;
    if (max < totptsvdim) max = totptsvdim;
  }
  if (L->work == NULL)  L->work = (double*) MALLOC(max * sizeof(double));
    
  PROTECT(all_res =  allocVector(VECSXP, sets));
  for (int set = 0; set < sets; set++){
    cov->base->set = set;
    if ((matrix = NCOL_OUT_OF(datasets) > 1)) break;
  }

  for (int set = 0; set < sets; set++){
    cov->base->set = set;
   int 
      ncol = NCOL_OUT_OF(datasets),
      nrow = NROW_OUT_OF(datasets);

    if (matrix) {
      PROTECT(res = allocMatrix(REALSXP, nrow, ncol));
    } else {
      PROTECT(res =allocVector(REALSXP, nrow));
    }

    get_logli_residuals(cov, L->work, REAL(res), modus);
    SET_VECTOR_ELT(all_res, set, res);
    UNPROTECT(1);
  }
  
  UNPROTECT(1);
  return all_res;
}

SEXP get_logli_residuals(SEXP model_reg) {
  START_GAUSS_INTERFACE;
  SEXP ans = get_logli_residuals(cov, LOGLI_RESIDUALS); // no PROTECT( needed
  return ans;
}

SEXP get_logli_wholetrend(SEXP model_reg) {
  //BUG;  printf("whats that good for?");
  START_GAUSS_INTERFACE;
  SEXP ans = get_logli_residuals(cov, LOGLI_WHOLETREND);// no PROTECT( needed
  return ans;
}



#define START_PREDICT							\
  model *conditioning = cov->key != NULL				\
    ? cov->key : cov->sub[PREDICT_CONDITIONING],			\
    *predict = cov->sub[PREDICT_PREDICT];				\
  assert(COVNR == PREDICT_CALL);					\
  assert(isnowProcess(conditioning) && conditioning->Ssolve != NULL &&	\
	 MODELNR(conditioning) == GAUSSPROC);				\
  assert(conditioning->Sfctn != NULL && conditioning->Slikelihood != NULL && \
	 isnowVariogram(conditioning->key == NULL			\
			? conditioning->sub[0] : conditioning->key));	\
  if (isnowProcess(predict))						\
    predict = predict->key == NULL ? predict->sub[0] : predict->key;	\
  assert(isnowVariogram(predict));					\
  assert(predict->Sfctn != NULL);					\
  assert(conditioning -> Slikelihood != NULL);				\
  GETSTORAGE(L, conditioning, likelihood);				\
  assert(L != NULL);							\
  listoftype *datasets = L->datasets;					\
  assert(L->datasets != NULL);						\
  int									\
  err = NOERROR,							\
    vdim = VDIM0,							\
    ncol = NCOL_OUT_OF(datasets),					\
    repet = ncol / vdim;

 

void get_fx(model *cov, double *v, int set) {
  // do NOT set cov->base->set
  int store = cov->base->set;
  cov->base->set = set;

  START_PREDICT;
  int pred_tot = Loctotalpoints(cov);					
  Long pred_totvdim = (Long) pred_tot * vdim;

  
  // currently unsused! kriging.variance?!
  // TO DO WITH 
  //
  // muss angenommen werden, dass alles gesetzt ist.
  assert (set >= 0 && set < LocSets(predict));
  int 
    betatot = L->cum_n_betas[L->fixedtrends],    
    pred_totvdimrepet = pred_tot * ncol;

  if (L->betas_separate) repet = 1;


 assert(false); // pred_tot f.s. falshc
  double *X = NULL;
  for (int k=0; k<pred_totvdimrepet; v[k++] = 0.0);
  if (L->ignore_trend) goto ErrorHandling; // exist
  
  
  if (!L->betas_separate && repet > 1) GERR("BUG");

  assert(false); // X -> TALLOC
  if ((X = (double*) MALLOC(pred_totvdim * sizeof(double))) == NULL) {
    // printf("%d %d \n", pred_tot, vdim); APMI(cov);
    err = ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  
  int r,m;
  for (int i=0; i<L->dettrends; i++) {
    assert(false); // naechstes Zeile sicher falsch
    FctnIntern(predict, L->det_effect[i],  L->det_effect[i], true, X);
    for (r = m = 0; r < repet; r++) {
      for (int k=0; k<pred_totvdim; k++) v[m++] += X[k];
    }
  }
  /*if(L->dettrends){
    for (z=r=0; r<repet; r++)
    for (int k=0; k<pred_totvdim; k++)
        v[z++] += L->YhatWithoutNA[cov->base->set][k];
  }*/
  for (int i=0; i<L->fixedtrends; i++) {
    assert(false);// naechstes Zeile sicher falsch
    FctnIntern(predict, L->fixed_effect[i],  L->fixed_effect[i], true, X);
    if (L->cum_n_betas[i+1] - L->cum_n_betas[i] != 1) BUG;
    double *betavec = L->betavec + L->cum_n_betas[i];
    for (r = m = 0; r < repet; r++) {
      double beta = *betavec;
      	//printf("beta = %10g %10g\n", beta, X[0]); BUG;
      for (int k=0; k<pred_totvdim; k++) v[m++] += X[k] * beta;
      if (L->betas_separate) betavec += betatot;
    }
  }
  
  ErrorHandling :
    cov->base->set = store;
  FREE(X);
  OnErrorStop(err, predict);
}


void get_F(model *cov, double *work, double *ans) {
  // do NOT set cov->base->set = 0;

   START_PREDICT;
  int pred_tot = Loctotalpoints(cov);					
Long pred_totvdim = (Long) pred_tot * vdim;

  assert(ans != NULL);
  
  int
    betas = L->cum_n_betas[L->fixedtrends],
    nrow = NROW_OUT_OF(datasets),
    //    ndata = nrow * ncol,
    totptsvdim = nrow * vdim; assert(false); // predict??!
  double *X = L->X[cov->base->set],
    *pres = ans,
    *data = SET_OUT_OF(datasets);
  
  //MEMCOPY(pres, data, ndata * sizeof(double));
  boxcox_trafo(P(GAUSS_BOXCOX), vdim, pres, nrow, repet);
  
  //   for (int i=0; i<ndata; i++) printf("%10g ", pres[i]);  BUG;
  //  printf("\n%d %d %d vdim=%d\n", L->ignore_trend, L->dettrends, L->fixedtrends, vdim);
                  
  if (L->ignore_trend) {
    return;
  }
                  
  bool delete_work;
  if ((delete_work = work == NULL)) {
    work = (double *) MALLOC(nrow * vdim * sizeof(double));
  }
  
  double *betavec = L->betavec;
  int r, z;
  
  if (L->dettrends) {
    
    for (int i=0; i<L->dettrends; i++) {      
      if (L->nas_det[i]) {
	assert(false); // stimmt nicht?!
	FctnIntern(cov, L->det_effect[i], L->det_effect[i], true,
		   work);
	for (z=r=0; r<repet; r++)
	  for (int k=0; k<totptsvdim; k++) pres[z++] += work[k];
      }
    }
    for (z=r=0; r<repet; r++)
      for (int k=0; k<totptsvdim; k++)
	pres[z++] += L->YhatWithoutNA[cov->base->set][k];
  }
  
  if (L->fixedtrends) {
    
    for (r=0; r<repet; r++, betavec += betas) {
      if (r==0 || L->betas_separate) {
	for (int i=0; i<totptsvdim; work[i++] = 0.0);
	//  for (int j=0; j<totptsvdim;j ++) printf("j=%d p %10g d %10g w %10g\n ", j, pres[j], data[j], work[j]);
	for (int i=0; i<betas; i++) {
	  double beta = betavec[i];
	  //	    printf("i=%d beta = %10g %d\n", i, beta, nrow);
	  for (int j=0; j<nrow; j++, X++) {	
	    //	      printf("i=%d %d %10g %10g\n", i, betas, beta, *X);
	    work[j] += *X * beta;
	  }
	}
      }
      
      for (int j=0; j<nrow; j++, data++, pres++) {
	*pres += work[j];
	//    printf("j=%d p %10g d %10g w %10g\n ", j, *pres, *data, work[j]);
      }
    }
  }
  
  // for (int j=0; j<nrow; j++, data++, pres++) printf("%10g %10g\n", *pres, *data);
  if (delete_work) FREE(work);
}




#define CHOLESKY_ONLY\
  solve_param Sparam;						\
  MEMCOPY(&Sparam, &(cov->base->global_utils.solve), sizeof(solve_param)); \
  Sparam.Methods[0] = Cholesky;						\
  Sparam.Methods[1] = NoFurtherInversionMethod

void gauss_trend(model *cov, double *v, int set, bool ignore_y) {
  //  printf("entering gauss trend\n");
  
  int store = cov->base->set;
  cov->base->set = set;

  START_PREDICT;
  int pred_tot = Loctotalpoints(cov);					
  Long pred_totvdim = (Long) pred_tot * vdim;

  assert (set >= 0 && set < LocSets(cov));
  int  betatot = L->cum_n_betas[L->fixedtrends];
  Long pred_totvdimrepet = (Long) pred_tot * ncol;
  
  for (Long k=0; k<pred_totvdimrepet; v[k++] = 0.0)
    ;

  // APMI(predict);
  TALLOC_XXX1(X, pred_totvdim);
  //  printf("Gauss trend %d %d %d bug:%d %d %d\n", pred_totvdim, L->ignore_trend, ignore_y,  !L->betas_separate && repet > 1, L->betas_separate, repet);
  
  if (L->ignore_trend) goto ErrorHandling; // exist
 
  int r,m;
  for (int i=0; i<L->dettrends; i++) {      
    FctnIntern(predict, L->det_effect[i], L->det_effect[i],
	       ignore_y, X);
    for (r = m = 0; r < repet; r++) {
      for (int k=0; k<pred_totvdim; k++) v[m++] += X[k];
    }
  }
  for (int i=0; i<L->fixedtrends; i++) {      
    FctnIntern(predict, L->fixed_effect[i],  L->fixed_effect[i],
	       ignore_y, X);
    if (L->cum_n_betas[i+1] - L->cum_n_betas[i] != 1) BUG;
    double *betavec = L->betavec + L->cum_n_betas[i];
    for (r = m = 0; r < repet; r++) {
      double beta = *betavec;
      //	printf("beta = %10g %10g\n", beta, X[0]); BUG;
      for (int k=0; k<pred_totvdim; k++) v[m++] += X[k] * beta;
      if (L->betas_separate) betavec += betatot;
    }
  }

 ErrorHandling :
  cov->base->set = store;

 END_TALLOC_XXX1;

// printf("end gauss trend %d\n", err);
  
  if (err != NOERROR) XERR(err);
}

/*
 if (global->krige.ret_variance) {
    //
    GERR("kriging variance cannot be returned, currently");
     // end !data_nas
  } else {



*/


void gauss_predict(model *cov, double *v) {
  // printf("entering gausss predict\n");
  // do NOT set cov->base->set
  globalparam *global = &(cov->base->global);
  assert(!global->krige.ret_variance);
  if (global->general.vdim_close_together) {
    ERR("'vdim_close_together' must be false for kriging and conditional simulation");
  }

  // printf("A\n");
  
  START_PREDICT;
  int pred_tot = Loctotalpoints(cov);					
  Long pred_totvdim = (Long) pred_tot * vdim;
 // printf("BA\n");

  int 
    Exterr = NOERROR,
    sets = LocSets(cov);
  if (sets != 1) ERR("only one data set allowed when predicting");
  assert(conditioning->Ssolve != NULL);
  
  int 
    spatialdim = Locspatialdim(cov),
    timespacedim = Loctsdim(cov),
    totptsY = LoctotalpointsY(cov), // conditioning points
    totptsYvdim  = totptsY * vdim,
    totptsYvdimSq  = totptsYvdim * vdim,
    totYvdimrepet = totptsY * vdim * repet;
  bool
    has_nas = L->data_nas[cov->base->set];

  CHOLESKY_ONLY;
 // printf("CBA\n");
  CovarianceMatrix(conditioning, false, L->C);  // tot x tot x vdim x vdim //sub -> cov 30.12.15
  // printf("BCDA\n");

  bool with_trend = true;
  if (with_trend) gauss_trend(cov, v, cov->base->set, true);//vor TALLOC sein!
  else for(int i=0; i<pred_totvdim; v[i++] = 0.0);
  // printf("BCDAX\n");


  //PMI0(cov);  printf("totYvdimrepet=%d %d %d, vdim=%d set=%d\n", totYvdimrepet, repet, ncol, vdim,cov->base->set);
  TALLOC_X1(residuals, totYvdimrepet); // X1 egal, da eh zu wenig

   // printf("BCDA\n");
//// printf("RESIFUALS %d\n", totYvdimrepet);

  
  get_logli_residuals(conditioning, NULL, residuals, LOGLI_RESIDUALS);
  //for (int i=0; i<201; i++) printf("%f ", residuals[i]); printf("\n");


  TALLOC_XXX1(cross, totptsYvdimSq);
  //  printf(" ACROSS tot = %d %ld\n", totptsYvdimSq, cross);
  TALLOC_XXX2(dummy, has_nas ? totptsYvdimSq : 1);

 // printf("WBA\n");
  double *ResiWithoutNA = NULL,
    *Ccur, *ptcross;
  int atonce, endfor,
    notnas = NA_INTEGER;
  if (has_nas) {
    atonce = 1;
    Ccur = L->Cwork;
    ptcross = dummy;
    ResiWithoutNA = L->Xwork;
   } else {
    notnas = totptsYvdim;
    atonce = repet;
    Ccur = L->C;
    ptcross = cross;
  }
  endfor = repet / atonce;
 
 // printf("BRA\n");

  for (int r=0; r<endfor; r++) { // over alle repetitions (in case of splitting)
    double *resi = residuals + r * totptsYvdim * atonce,
      *pv0 = v + r * atonce * pred_totvdim;
    if (has_nas) {
      assert(ResiWithoutNA != NULL);
      notnas = matrixcopyNA(ResiWithoutNA, resi, resi, totptsYvdim, 1, 0);

      //  for (int i=0; i<totptsYvdim*totptsYvdim ; i++) {
      ///	if (i % totptsYvdim ==0)printf("\n");
      //	printf("%d>%e ", i, L->C[i]);
      // } printf("\nEnde L->C \n\n");
         
      SqMatrixcopyNA(Ccur, L->C, resi, totptsYvdim);
    } else ResiWithoutNA = resi;

    
    //printf("hasnas = %d %d %d %ld\n", has_nas, notnas, totptsYvdim, ResiWithoutNA);

    //   printf("%d\n", resi - residuals);
    //   for (int i=0; i<notnas; i++) printf("%d>%5.4f ", i, residuals[i]);
    // printf("\nende Resi\n");
  
 
    //  for (int i=0; i<notnas*notnas ; i++) {
    //	if (i % notnas ==0)printf("\n");
	//	printf("%5.2f ", Ccur[i]);
    //    }
    //  printf("\nende CCur\n%d notnas=%d\n", totptsYvdim, notnas);
  
    
    //  for (int i=0; i<notnas; i++) printf("%d>%5.4f ", i, ResiWithoutNA[i]);
    //  printf("\nende Resi\n");

    
    assert(conditioning->Ssolve != NULL);
    Exterr = Ext_solvePosDefSp(Ccur, notnas, 
			       true, // except for intrinsic Kriging
			       //  as if ordinary Kriging (p 167; 265 in Chiles
			       ResiWithoutNA, atonce, NULL,
			       conditioning->Ssolve,
			       &Sparam);
    if (Exterr != NOERROR)
      GERR2("In RandomFieldsUtils: %.200s (error = %d)", cov->Ssolve->err_msg,
	    Exterr);

    //   for (int i=0; i< notnas; i++) printf("%i:%5.2f ", i, ResiWithoutNA[i]); printf("\n");

      
    for (int i_row=0; i_row<pred_tot; i_row++, pv0++) { // over all locations
     //for ease: simple kriging considered:
      // calculated above: Z = (\Sigma_{y_jy_j})J^{1} %*% residuals
      // to be calulated S_i:=((\Sigma_{x_i, y_j}))_j,
      // which is a lying vector.
      // Technically faster is to use is S_i^T, so thats why Covariance*T.
      // In definition of CovVario: y is recycled. Here: x must be recycled.
      // So, CovarianceT also swaps role of x and y

       
      CovarianceT(predict, i_row, cross);
     
      if (has_nas) matrixcopyNA(dummy, cross, resi, totptsYvdim, vdim, 0);

      //printf("ptcross = %ld %ld\n", ptcross, cross);

      // printf("\ni=%d\n", i_row); for (int i=0; i< notnas; i++) printf("%d/%e ", i, cross[i]); printf("\n");

      double *pv = pv0;
      double *rWNA = ResiWithoutNA;
      for (int s=0; s<atonce; s++, rWNA += notnas) { // over indep. realisations
	double *ptc = ptcross;
	//warum nicht ResiWithoutNA statt residuak? -> geaendert 30.12.20
	//double factor = residuals[totptsYvdim * s + totptsY * j+i_row];

	for (int j=0; j<vdim; j++, ptc += notnas, pv += pred_tot) {//multivar.
	  // influence of the different variates of position i

	  // printf("\nbeor vdm=%d, p_tot=%d nas=%d notnas=%d=%d  i=%d j=%d s=%d r=%d %e\n", vdim, pred_tot, L->nas_fixed[cov->base->set], notnas, totptsYvdim- L->nas_fixed[cov->base->set], i_row, j, s,r, *pv);
	  //double sum = 0.0; for (int i=0; i<notnas; i++) sum += ptcross[i] * ResiWithoutNA[i];

	  //	  printf("r=%d i=%d , s=%d, j=%d, ptcros %ld , resi %ld, notnas %d\n", r, i_row, s, j, ptc, rWNA, notnas);

	  (*pv) += SCALAR_P(ptc, rWNA, notnas);

	  //	  printf("i=%d j=%d s=%d %d %e \n", i_row, j, s, r, *pv);
	}// j, vdim
      } // s, atonce
    } // i_row
  } // r, repeated, not atonve 
  END_TALLOC_XXX1;
  END_TALLOC_XXX2;
  END_TALLOC_X1;

  //  printf("done\n");
  
 ErrorHandling :
  cov->base->set = 0;

  if (err != NOERROR) XERR(err);
}

//////////////////////////////////////////////////////////////////////
// Gauss process

SEXP simple_residuals(SEXP model_reg){
  START_GAUSS_INTERFACE;
  likelihood_facts *facts = &(L->facts);
  assert(cov->Ssolve != NULL);
  
  assert(!L->ignore_trend);
  
  listoftype *datasets = L->datasets;
  if (facts->globalvariance && facts->pt_variance != NULL) 
    *(facts->pt_variance) = 1.0;
 
  int i, j, 
    fx_notnas = 0, // number of betas that do not have nas
    sets = LocSets(calling),  
    err = NOERROR,
    Exterr = NOERROR,
    vdim = VDIM0,
    betatot = L->cum_n_betas[L->fixedtrends],
    betaSq = betatot * betatot;
  //if (VDIM1 != 1) BUG;

  CHOLESKY_ONLY;
  for (i=0; i<L->fixedtrends; i++) {
    fx_notnas += (!L->nas_fixed[i]) * (L->cum_n_betas[i+1] - L->cum_n_betas[i]);
    // printf("fx i=%d %d %d %d\n", i, L->nas_fixed[i], L->cum_n_betas[i+1], L->cum_n_betas[i]);
  }
  for (i=0; i<betaSq; L->XtX[i++] = 0.0);
  for (i=0; i<betatot; L->XCY[i++] = 0.0);
  for (int set=0; set < sets; set++) {
    cov->base->set = set;
    //  Loctotalpoints, Loc etc haengen alle von set ab !
    int k, m, 
      ncol =  NCOL_OUT_OF(datasets),
      repet = ncol / vdim,
      nrow = NROW_OUT_OF(datasets),
      ndata = nrow * ncol,
      totpts = LoctotalpointsY(calling),
      totptsvdim = totpts * vdim;
    assert(nrow == totpts); 
    double 
      *Xdata = L->X[set] + betatot * totptsvdim;
  
    if (L->dettrend_has_nas || L->nas_boxcox) {
      MEMCOPY(Xdata, SET_OUT_OF(datasets), ndata * sizeof(double));
      if (L->nas_boxcox) 
	boxcox_trafo(P(GAUSS_BOXCOX), vdim, Xdata, nrow, repet);
    }
   
    for (i=0; i<L->fixedtrends; i++) {
      if (L->nas_fixed[i]) {
	double *Lx = L->X[set] + L->cum_n_betas[i] * totptsvdim,
	  end = totptsvdim * (L->cum_n_betas[i+1] - L->cum_n_betas[i]);
	for (m=0; m<end; Lx[m++] = 0.0);
      }
    }

    if (L->random > 0) BUG; // to do
 
    double *Xcur,
      *dataWithoutNA = NULL;
    int atonce, endfor,
      notnas = NA_INTEGER;
    if (L->data_nas[set]) {
      atonce = 1;
      Xcur = L->Xwork;
    } else {
      atonce = repet;
      Xcur = L->X[set];
      dataWithoutNA = Xdata;
      notnas = totptsvdim;
    }
    endfor = repet / atonce;
    
 
    for (int r=0; r<endfor; r++) {
      if (L->data_nas[set]) {
	double *datacur = Xdata + r * totptsvdim;
	notnas = matrixcopyNA(Xcur, L->X[set], 
			      datacur, totptsvdim, betatot, atonce);
	dataWithoutNA = Xcur + notnas * betatot;
      }

      if (L->fixedtrends) {
	// gleoch alles mit allem Multiplizieren, um's einfach zumachen.
	// relevantes wird dann aus L->XitXi rausgeholt.
	matmulttransposed(Xcur, Xcur, L->XitXi, notnas, betatot, betatot);
	// bei simple residuals muessen alle spalten und Zeilen
	// mit nas_fixed raus:
	for (int i1=k=0; i1<L->fixedtrends; i1++) { 
	  if (L->nas_fixed[i1] == 0) {
	    int end_i = L->cum_n_betas[i1 + 1];
	    for (int i2=L->cum_n_betas[i1]; i2<end_i; i2++) {
	      for (int j1=0; j1<L->fixedtrends; j1++) { 
		if (L->nas_fixed[j1] == 0) {
		  int end_j = L->cum_n_betas[j1 + 1];
		  for (int j2=L->cum_n_betas[j1]; j2<end_j; j2++) {
		    int idx = i2 * betatot + j2;
		    L->XtX[k++] += atonce * L->XitXi[idx];
		  } // j2
		}
	      } // j1
	    } // i2
	  }
	} // i1
	
	MEMCOPY(L->sumY, dataWithoutNA, sizeof(double) * notnas);
	k = notnas;
	for (int rr=1; rr<atonce; rr++) {
	  for(int d=0; d<notnas; d++) {
	    L->sumY[d] += dataWithoutNA[k++];
	  }
	}
	
	double *dummy = L->betavec;
	matmult(L->sumY, Xcur, dummy, 1, notnas, betatot);
	for (i=0; i<betatot; i++) L->XCY[i] += dummy[i];
      }
    } // for ... end !data_nas
  } // set


  // PMI(cov);
  
  if (L->fixedtrends) {
    // XCY, etc NUR FUER FIXED TREND
    double *beta = L->betavec; // !!!
    for (i = 0; i<betatot;  beta[i++] = 0.0);

    if (fx_notnas > 0) {
      int k=0;
      for (int j1=0; j1 < L->fixedtrends; j1++) { 
	if (L->nas_fixed[j1] == 0) {
	  int end=L->cum_n_betas[j1 + 1];
	  for (int j2=L->cum_n_betas[j1]; j2<end ; j2++) {
	    beta[k++] = L->XCY[j2];
	  } // j2
	}
      }
      //

      Exterr = Ext_solvePosDefSp(L->XtX, fx_notnas, true, beta, 1, NULL,  
				 cov->Ssolve, &Sparam);
      if (Exterr != NOERROR)
	GERR2("In RandomFieldsUtils: %.200s (error = %d)", cov->Ssolve->err_msg,
	      Exterr);


      
      if (fx_notnas < betatot) {
	int b = L->fixedtrends - 1;
	i = betatot - 1; 
	j = fx_notnas - 1;
	//	printf("i=%d j=%d b=%d\n", i, j, b);
	while (b>=0 && i >= 0) {	 
	  int bi = L->cum_n_betas[b + 1] - L->cum_n_betas[b];
	  //	  printf("bi=%d\n", bi);
	  if (L->nas_fixed[b]) {
	    for (k=0; k<bi; k++) beta[i--] = 0.0;
	  } else for (k=0; k<bi; k++) {
	      //     printf("k=%d i=%d j=%d\n", k, i, j);
	      assert(i >= 0 && j >= 0);
	      beta[i--] = beta[j--];
	    }
	  b--;
	  //  printf("i=%d j=%d\n", i, j);
	  assert(i>=j);
	}
        assert (i == -1 && j == -1 && b == -1);
      } else assert(fx_notnas == betatot);
    }
  }
  
 ErrorHandling:
  cov->base->set = 0;
  if (err != NOERROR) XERR(err);

  return get_logli_residuals(model_reg);
}


//////////////////////////////////////////////////////////////////////
// Gauss process

// likelihood / density
void gaussprocessDlog(double VARIABLE_IS_NOT_USED *x, INFO, model *cov, 
		      double *v){
  cov->base->set = 0;
  assert(isnowProcess(cov));
  model *calling = cov->calling;
  //*sub = cov->key == NULL ? cov->sub[0] : cov->key;  
  assert(calling != NULL && DefList[CALLINGNR].cov == likelihood);
getStorage(L ,   likelihood);
  assert(L != NULL);
  listoftype *datasets = L->datasets;
  assert(!L->ignore_trend);
  likelihood_facts *facts = &(L->facts);
  if (facts->globalvariance && facts->pt_variance != NULL) 
    *(facts->pt_variance) = 1.0;
  double 
    proj = 0.0,
    YCinvY = 0.0,
    logdettot = 0.0;  
  int i, j, repet,
    sets = LocSets(calling),  
    err = NOERROR,
    Exterr = NOERROR,
    variables = 0,
    vdim = VDIM0,
    betatot = L->cum_n_betas[L->fixedtrends],    
    betaSq = betatot * betatot;
  // PMI(cov);
  //  if (VDIM1 != 1) BUG;
  CHOLESKY_ONLY;

  *v = 0.0;
  for (i=0; i<betaSq; L->XtX[i++] = 0.0);
  for (i=0; i<betatot; L->XCY[i++] = 0.0);
  for (int set = 0; set < sets; set++) {
    cov->base->set = set;
    //  Loctotalpoints, Loc etc haengen alle von set ab !    
    int k, m, p,
      ncol =  NCOL_OUT_OF(datasets);    
    assert(!L->betas_separate || sets == 1);
    repet = ncol / vdim;
    int nrow = NROW_OUT_OF(datasets), 
      ndata = nrow * ncol,
      totpts = LoctotalpointsY(calling),
      totptsvdim = totpts * vdim;
    double 
      logdet,
      *val = L->C,
      *YhatWithoutNA = L->YhatWithoutNA[set],
      *Xdata = L->X[set] + betatot * totptsvdim;
  
    if (L->dettrend_has_nas || L->nas_boxcox) {
      MEMCOPY(Xdata, SET_OUT_OF(datasets), ndata * sizeof(double));
 
      if (L->nas_boxcox) {
	for (p = m = 0; m < repet; m++) {
	  for (k = 0; k < totptsvdim; k++) {
	    Xdata[p++] -= YhatWithoutNA[k];
	  }
	}
      }
    
      for (i=0; i<L->dettrends; i++) {
	if (L->nas_det[i]) {
	  FctnIntern(cov, L->det_effect[i], L->det_effect[i], false, val);  
	  for (p = m = 0; m < repet; m++) {
	    for (k = 0; k < totptsvdim; k++) Xdata[p++] -= val[k];
	  }
	}
      }

      if (L->nas_boxcox) 
	boxcox_trafo(P(GAUSS_BOXCOX), vdim, Xdata, nrow, repet);

    } // (L->dettrend_has_nas || L->nas_boxcox)

    for (i=0; i<L->fixedtrends; i++) {
      //      printf("nas->fixed: i=%d %d\n", i, L->nas_fixed[i]);
      if (L->nas_fixed[i]) {
	FctnIntern(cov, L->fixed_effect[i], L->fixed_effect[i], false,
		   L->X[set] + L->cum_n_betas[i] * totptsvdim);
	// "Lueckentext" auffuelen !!
      }
    }

    if (L->random > 0) BUG; // to do
    if (facts->globalvariance && facts->pt_variance != NULL)
      *(facts->pt_variance) = 1.0;

    ///   PMI0(cov);
    
    CovarianceMatrix(cov, false, L->C); //sub -> cov 30.12.15
  
   
    //    printf("%10g %10g %10g %10g %10g\n", Xdata[0], Xdata[1], Xdata[2], Xdata[3], Xdata[4]); 		      

    double *Xcur, *Ccur, 
      *dataWithoutNA = NULL;
    int atonce, endfor, XYcols, 
      notnas = NA_INTEGER;
    if (L->data_nas[set]) {
      atonce = 1;
      Ccur = L->Cwork;
      Xcur = L->Xwork;
    } else {
      notnas = totptsvdim;
      atonce = repet;
      Ccur = L->C;
      Xcur = L->X[set];
      dataWithoutNA = Xdata;
    }
    endfor = repet / atonce;
    XYcols = betatot + atonce;

    // printf("%10g\n", L->X[set][0]); BUG; 

    for (int r=0; r<endfor; r++) {
      if (L->data_nas[set]) {
	double *datacur = Xdata + r * totptsvdim;
	notnas = matrixcopyNA(Xcur, L->X[set], 
			      datacur, totptsvdim, betatot, atonce);
	SqMatrixcopyNA(Ccur, L->C, datacur, totptsvdim);
	dataWithoutNA = Xcur + notnas * betatot;
      }
 
      variables += notnas * atonce;
      MEMCOPY(L->CinvXY, Xcur, notnas * XYcols * sizeof(double));
      assert(cov->Ssolve != NULL);

      //      BUG;
      // kosten vom aufruf : 2 / 1
      // biwm2 3 / 8
      //      printf("XXX= %d\n", 4);
      //
       
      Exterr = Ext_solvePosDefSp(Ccur, notnas, 
				 true, // except for negative definite function
				 L->CinvXY, XYcols, 
				 &logdet, cov->Ssolve, &Sparam);
      if (Exterr != NOERROR)
	GERR2("In RandomFieldsUtils: %.200s (error = %d)", cov->Ssolve->err_msg,
	      Exterr);

      
      logdettot += logdet * atonce;
      double *CinvY = L->CinvXY + betatot * notnas;
      int end_i = notnas * atonce;

      //printf("\nCinvXY");for (int i=0; i<5; i++) printf("%f ", L->CinvXY[i]);
      
      // for (i=0; i<end_i; i++) YCinvY += CinvY[i] * dataWithoutNA[i];
      YCinvY += SCALAR_P(CinvY, dataWithoutNA, end_i);

      if (L->fixedtrends) {
	if (L->betas_separate) {
	  assert(sets == 1);
	  matmulttransposed(Xcur, L->CinvXY, L->XtX, notnas, betatot, betatot);
	  matmulttransposed(L->CinvXY, dataWithoutNA, L->XCY, notnas,
			    betatot,atonce);
	} else {
	  //BUG;
	  matmulttransposed(Xcur, L->CinvXY, L->XitXi, notnas, betatot,betatot);
	  for (i=0; i<betaSq; i++) L->XtX[i] += atonce * L->XitXi[i];

	  // printf("atonce=%d %10g xcur=%10g %10g\n", atonce, L->XitXi[0], Xcur[0], L->CinvXY[0]);
	  
	  MEMCOPY(L->sumY, dataWithoutNA, sizeof(double) * notnas);

	  int kk = notnas;
	  for (int rr=1; rr<atonce; rr++) {
	    for(int d=0; d<notnas; d++) {
	      L->sumY[d] += dataWithoutNA[kk++];
	    }
	  }
	  
	  double *dummy = L->betavec;
	  matmult(L->sumY, L->CinvXY, dummy, 1, notnas, betatot);	  
	  for (i=0; i<betatot; i++) L->XCY[i] += dummy[i];
	}
      } // fixedtrends	
    } // r, repet	  
  } // set
  
  if (L->fixedtrends) {
    // XCY, etc NUR FUER FIXED TREND
    assert(!L->betas_separate || sets == 1);
    double *beta = L->betavec; // !!!
    int 
      all_betatot = betatot,
      ncol =  NCOL_OUT_OF(datasets);
    repet = ncol / vdim;
  
    if (L->betas_separate) {
      all_betatot *= repet;
    }
  
   
    MEMCOPY(beta, L->XCY, sizeof(double) * all_betatot);
    assert(!L->betas_separate || sets == 1);
    assert(cov->Ssolve != NULL);

    //    for (int h=0; h<betatot * betatot; h++) printf("%10g ", L->XtX[h]);
    //    printf("\n>> %d %d %10g\n",  betatot, L->betas_separate, beta[0]);

    Exterr = Ext_solvePosDefSp(L->XtX, betatot, true, beta, 
			       L->betas_separate ? repet : 1,
			       NULL, 
			       cov->Ssolve, &Sparam);
    if (Exterr != NOERROR)
      GERR2("In RandomFieldsUtils: %.200s (error = %d)", cov->Ssolve->err_msg,
	    Exterr);
    
    //   printf("%.50s AC %10g %10g %10g %10g %d done\n", NAME(cov), beta[0], beta[1], beta[2], beta[3], all_betatot); BUG;
    for (j=0; j<all_betatot; j++) {
      proj += L->XCY[j] * beta[j];
    }      
    MEMCOPY(v + 1 + facts->globalvariance, beta, all_betatot * sizeof(double));

    if (L->betas_separate) for(i=0; i<betatot; i++) *(L->where_fixed[i]) =RF_NA;
    else for (i=0; i<betatot; i++) *(L->where_fixed[i]) = beta[i];
  }

  *v = -0.5 * (logdettot + variables * LOG(TWOPI));

  //  printf(">> v =%10g %10g #=%d 2pi^=%10g sq=%10g %10g %d\n", *v, 0.5 *logdettot, variables , -0.5 * variables * LOG(TWOPI), YCinvY, proj, facts->globalvariance); 
 
  double delta;
  delta = YCinvY - proj;
  if (facts->globalvariance) {
    //   printf("%10g delta=%10g %d %10g \n", *v, delta, variables, 0.5 * variables * LOG(delta)); BUG;
    
    v[1] = delta / variables;
    *v -= 0.5 * variables * (1 + LOG(v[1]));
    if (facts->pt_variance != NULL) *(facts->pt_variance) = v[1];
  } else {
    //    printf("delta = %10g, %10g YCY=%10g proj=%10g\n", delta, -0.5 * delta, YCinvY,proj);
    *v -= 0.5 * delta;
  }
  //  printf("*v=%10g %10g %d\n", *v, P(GAUSS_BOXCOX)[0], facts->globalvariance);
  
  if (R_FINITE(P(GAUSS_BOXCOX)[0])) {
    //printf("*v=%10g %10g %d\n", *v, P(GAUSS_BOXCOX)[0], facts->globalvariance);BUG;

    for (j=0; j<vdim; j++) {
      double lambda  = P(GAUSS_BOXCOX)[2 * j];
      if (lambda != 1.0) {
	double 
	  lambdaM1 = lambda - 1.0,
	  mu = P(GAUSS_BOXCOX)[2 * j + 1],
	  sum = 0.0;
	for (int set = 0; set < sets; set++) {
	  cov->base->set = set;
	  int 
	    ncol =  NCOL_OUT_OF(datasets),   
	    nrow = NROW_OUT_OF(datasets), 
	    totptsvdim = nrow * vdim;
	  repet = ncol / vdim;
	  double 
	    *data = SET_OUT_OF(datasets) + nrow * j;
	  for (int r=0; r < repet; r++, data += totptsvdim) {
	    for (int ii=0; ii<nrow; ii++) {
	      double dummy = data[ii];
	      if (!ISNAN(dummy))
		sum += LOG(dummy + mu);
	    }
	  }
	}
	*v += lambdaM1 * sum;
      }
    }
  }
  //  printf("      *v=%10g %10g\n", *v, P(GAUSS_BOXCOX)[0]);

  
  //  APMI(cov);  BUG;

 ErrorHandling:
  cov->base->set = 0;

  if (err != NOERROR) XERR(err);
   
}

// v ok
void PutGlblVar(int *reg, double *var) {
  assert (*reg >= 0 && *reg <= MODEL_MAX);
  model *cov = KEY()[*reg];
  assert(cov != NULL && equalsnowInterface(cov));
  model *sub = cov->key == NULL ? cov->sub[0] : cov->key; 
  likelihood_storage *L;
  if (sub == NULL || !isnowProcess(sub) || (L = sub->Slikelihood) == NULL) BUG;
  likelihood_facts *facts = &(L->facts);
  if (facts->pt_variance != NULL) {
    // == NULL may happen if globalvariance is true by user, but not
    // reflected in the model
    *(facts->pt_variance) = *var;
  }
}

SEXP get_likeliinfo(SEXP model_reg) {
 START_GAUSS_INTERFACE;
  likelihood_facts *facts = &(L->facts);
#define nn 5
  const char *names[nn] = 
    {"betas", "betanames", 
     "estimate_variance", "sum_not_isna_data", 
     "betas_separate"};
  listoftype *datasets = L->datasets;
  int k, 
    sets = LocSets(cov),  
    sum_not_isna_data = 0,
    betas = L->cum_n_betas[L->fixedtrends];
  SEXP namevec, ans, betanames;

  for (int set = 0; set < sets; set++){
    cov->base->set = set;
    int
      ncol = NCOL_OUT_OF(datasets),
      nrow = NROW_OUT_OF(datasets),
      ndata = nrow * ncol;
    sum_not_isna_data += ndata - L->data_nas[set];
  }
  
  PROTECT(ans = allocVector(VECSXP, nn));
  PROTECT(namevec = allocVector(STRSXP, nn));
  for (k=0; k<nn; k++) SET_STRING_ELT(namevec, k, mkChar(names[k]));

  //  printf("betas %d\n", betas);
  PROTECT(betanames = allocVector(STRSXP, betas));
  for (k=0; k<betas; k++) {
    //printf("%d %.50s<\n", k, L->betanames[k]);
    SET_STRING_ELT(betanames, k, mkChar(L->betanames[k]));
  }

  k = 0;  
  SET_VECTOR_ELT(ans, k++, ScalarReal(betas));
  SET_VECTOR_ELT(ans, k++, betanames);
  SET_VECTOR_ELT(ans, k++, ScalarLogical(facts->globalvariance));
  SET_VECTOR_ELT(ans, k++, ScalarInteger(sum_not_isna_data));
  SET_VECTOR_ELT(ans, k++, ScalarLogical(L->betas_separate));
  assert(k == nn);


  setAttrib(ans, R_NamesSymbol, namevec);
  UNPROTECT(3);
 
  cov->base->set = 0;
  return ans;
}


int countbetas(model *cov, double ***where) {
  defn *C = DefList + COVNR;
  int i,
    sum = 0,
    kappas = C->kappas;

  for(i=0; i<kappas; i++) {
    if (cov->kappasub[i] != NULL) continue;
    if (equalsnowTrendParam(cov, i)) {
      assert(C->kappatype[i] == REALSXP);
      if (!PisNULL(i)) {
	int total = cov->ncol[i] * cov->nrow[i];
	double *p = P(i);
	if (ISNAN(p[0])) {
	  sum += total;
	  for (int j=0; j<total; j++) { // 0 for 'where'
	    if (!ISNAN(p[j])) ERR("trend parameters must be all NA or none");
	    if (where != NULL) {
	      **where = p + j;
	      (*where)++;	     
	    }	    
	  }
	} else {
	  for (int j=1; j<total; j++) {
	    if (ISNAN(p[j])) ERR("trend parameters must be all NA or none");
	  } // j
	} // not (ISNAN(p[0])) 
      } // (!PisNULL(i)) 
    } // is Trend
  } // for kappa
  //  printf("sum = %d\n", sum);
  return sum;
}



void AbbrBeta(char *Old, char *abbr) {
  int len = 6;// global->fit.lengthshortname / 3,
  if (Old == NULL) {
    STRNCPY(abbr, "unk.", len);
    return;
  }

    char *old = Old;
  if (old[0] == '.') old++;

  // printf("%d\n", global->fit.lengthshortname ); BUG;
  
  int nold = STRLEN(old),
    nabbr = len - 1;

  if (nold <= len) {
    STRNCPY(abbr, old, len);
    return;
  }
  abbr[0] = old[0];
  while (nabbr >= 1 && nabbr < nold) { 
    char b = old[nold];
    if (b=='a' || b=='A' || b=='e' || b=='E' || b=='i' || b=='I' ||
        b =='o' || b=='O' || b=='u' || b=='U') nold--;
    else abbr[nabbr--] = old[nold--];
  }
  if (nabbr > 1) {
    assert(nabbr==0 || nold == nabbr);
    for (int i=2; i<=nold; i++) abbr[i] = old[i];
  }
  abbr[len] = '\0';
}
  
  
void GetBeta(model *cov, likelihood_storage *L, int *neffect)  {
  globalparam *global = &(cov->base->global);
  if (isnowProcess(cov)) {
    int nas = ((bool) ISNA((P(GAUSS_BOXCOX)[0]))) +
      ((bool) ISNA((P(GAUSS_BOXCOX)[1])));
    assert(cov->key == NULL);
    assert(!PisNULL(GAUSS_BOXCOX));
    if (nas > 0) (*neffect)++;
    GetBeta(cov->sub[0], L, neffect);
    return;
  }

  int i,
    n = COVNR == PLUS ? cov->nsub : 1;
  if (*neffect >= MAX_LIN_COMP) ERR("too many linear components");
  bool 
    plus = COVNR == PLUS;
  likelihood_facts *facts = &(L->facts);

  //  printf("entering getbeta\n");  pmi(cov, 0);

  for (i=0; i<n; i++) {
    model *component = plus ? cov->sub[i] : cov;
    //   pmi(component, 0);
    if (MODELNR(component) == PLUS) {
      //   printf("recursive\n");
      GetBeta(component, L, neffect);
      //  printf("end recursive\n");
      continue;
    }

    // printf("%.50s %d %ld\n", NAME(cov), facts->effect[*neffect], (Long) L);
    
    //  if (L->effect[*neffect] > LastMixedEffect) { } else 

    //printf("neffect %d %d\n", *neffect, facts->effect[*neffect]);
    if (facts->effect[*neffect]  == DetTrendEffect) {     
      L->det_effect[L->dettrends++] = component;
    } else if (facts->effect[*neffect] == FixedEffect) {
      int b = UNSET;
      L->cum_n_betas[L->fixedtrends + 1] = L->cum_n_betas[L->fixedtrends];
      L->fixed_effect[L->fixedtrends++] = component;
      if (MODELNR(component) == MULT) {
	for (int j=0; j<component->nsub; j++) {
	  if ((b = countbetas(component->sub[j], NULL)) > 0) {
	    break;
	  }
	}
      } else b = countbetas(component, NULL);	 

      if (b <= 0) {
	(*neffect)++;
	continue;
      }
      //ERR("fixed effect expected, but corresponding NA not found");
      int base =  L->cum_n_betas[L->fixedtrends];
      L->cum_n_betas[L->fixedtrends] += b;  
      if (L->maxbeta < b) L->maxbeta = b;

      model *comp = component;
      int ii = -1,
	nr = MODELNR(component);
      if (nr == MULT) {
	for (ii=0; ii<comp->nsub; ii++) {
	  model *subcomp = comp->sub[ii];
	  if (MODELNR(subcomp) == CONST && ISNA(PARAM0(subcomp, CONST_C))){
	    comp = comp->sub[ii==0 && comp->nsub >= 2];
	    break;
	  }
	}
      }
      if (isDollar(comp)) comp = comp->sub[0]; // jump just for a nice name
      char abbr[LENMSG];
      int bytes = (1 + global->fit.lengthshortname) * sizeof(char);
      //      PMI0(comp);
      //printf("%d %.50s %.50s %d\n", *neffect, NAME(component), NAME(comp),2);
      //      printf("%d\n",  PARAMisNULL(comp, COVARIATE_NAME));
      if (nr == COVARIATE)
	AbbrBeta(PARAMisNULL(comp, COVARIATE_NAME) ? NULL
		 : PARAM0CHAR(comp, COVARIATE_NAME), abbr);

      else if (nr == CONST) SPRINTF(abbr, "const");

      else if (nr == SHAPE_FCT && comp->kappasub[SHAPE_FCT_MEAN] != NULL) 
	AbbrBeta(NICK(comp->kappasub[SHAPE_FCT_MEAN]) + 2, abbr);

      else if (nr == MULT) {
	if (MODELNR(comp) == PROJ_MODEL) {
	  if (PARAMisNULL(comp, PROJ_NAME))
	    SPRINTF(abbr, "proj%d", PARAM0INT(comp, PROJ_PROJ));
	  else AbbrBeta(PARAMisNULL(comp, PROJ_NAME) ? NULL
			: PARAM0CHAR(comp, PROJ_NAME), abbr);
	}
	else if (MODELNR(comp) == SHAPE_FCT &&
		 comp->kappasub[SHAPE_FCT_MEAN] != NULL)
	  AbbrBeta(NICK(comp->kappasub[SHAPE_FCT_MEAN]) + 2, abbr);
	else AbbrBeta(NICK(comp) + 2, abbr);
      }
      
      else AbbrBeta(NICK(comp) + 2, abbr);
      //      printf("%.50s<\n", abbr);
      //  if (L->betanames[*neffect] != NULL) {printf("%.50s %d %.50s\n", NAME(cov),  *neffect, L->betanames[*neffect]); crash();} else  printf("%.50s %d NULL\n", NAME(cov),  *neffect); 
      
      assert(L->betanames[base] == NULL);
      if (b == 1) {
	L->betanames[base] = (char*) MALLOC(bytes);
	// printf("base = %s\n", abbr);
	SPRINTF(L->betanames[base], "%.50s", abbr);
      } else {
	for (int jj=0; jj<b; jj++) {
	  L->betanames[base + jj] = (char*) MALLOC(bytes);
	  SPRINTF(L->betanames[base + jj], "%.50s.%d", abbr, jj);
	}
      }

      //printf("----> %d %.50s\n", base, L->betanames[base]);

    } 
    (*neffect)++;
  }

  //printf("%d\n", L->cum_n_betas[L->fixedtrends]);
  //BUG;
  
}




void GetBeta(model *cov, likelihood_storage *L, int *neffect,
	     double ***where)  {
  if (isnowProcess(cov)) {
    int nas = ISNA((P(GAUSS_BOXCOX)[0])) + ISNA((P(GAUSS_BOXCOX)[1]));
    assert(cov->key == NULL);
    assert(!PisNULL(GAUSS_BOXCOX));
    if (nas > 0) (*neffect)++;
    GetBeta(cov->sub[0], L, neffect, where);
    return;
  }
  int i,
    n = COVNR == PLUS ? cov->nsub : 1;
  bool 
    plus = COVNR == PLUS;
  likelihood_facts *facts = &(L->facts);
  for (i=0; i<n; i++) {
    model *component = plus ? cov->sub[i] : cov;
    if (MODELNR(component) == PLUS) {
      GetBeta(component, L, neffect, where);
      continue;
    }
    if (facts->effect[*neffect] == FixedEffect) {
      if (MODELNR(component) == MULT) {
	for (int j=0; j<component->nsub; j++) {
	  if ((countbetas(component->sub[j], where)) > 0) {
	    break;
	  }
	}
      } else countbetas(component, where);	  
    }
    (*neffect)++;
  }
}


#define MAX_TOTALPTS 10000
int struct_gauss_logli(model *cov) {
  // APMI(cov);
  cov->base->set = 0;
  globalparam *global = &(cov->base->global);
  assert(isnowProcess(cov));
  model 
    *calling = cov->calling,
    *sub = cov->key != NULL ? cov->key : cov->sub[0];

  /*
    if (false) { // Beschleunigung, insb. im Genetik-Bereich
    // gleich ganz allgemein implementieren??
    if (isCartesian(SUB)) {
    if (isXonly(SUB)) {
    if (isIsotropic(SUB)) {
    if (Loctotalpoints(cov) < MAX_TOTALPTS) {
    // nur distances abspeichern
    // in ownloc
    }
    } else {
    if (Loctotalpoints(cov) <
    MAX_TOTALPTS / SQRT((double) dim)) {
    // distance vectors abspeochen
    }
    }
    }
    } else {
    NotProgrammedYet("non-cartesian systems")
    // to do 
    }
    }
  */

  int i, k,  // facts->neffect, 
    betatot,
    err = NOERROR,
    ne = 0,
    //n = SUBNR == PLUS ? sub->nsub : 1,
    sets = LocSets(cov),  
    maxrepet = 0,
    maxndata = 0,
    maxtotpts = 0,
    vdim = VDIM0,
    two_vdim = 2 * vdim,  
    n = vdim * vdim,
    dummy_fixed=0, 
    dummy_det = 0;
  usr_bool 
    global_var = CALLINGNR == LIKELIHOOD_CALL 
    ? (usr_bool) PARAM0INT(calling, LIKELIHOOD_NA_VAR)
    : False;
  Long max_total_data_Sq;
  
  assert(calling != NULL);
  if (CALLINGNR != LIKELIHOOD_CALL && CALLINGNR != LINEARPART_CALL &&
      CALLINGNR != PREDICT_CALL) BUG;
  NEW_STORAGE(likelihood);
getStorage(L ,   likelihood);
  likelihood_NULL(L);
  likelihood_facts *facts = &(L->facts);
  likelihood_facts_NULL(facts);
  L->ignore_trend = 
    CALLINGNR == LIKELIHOOD_CALL && PARAM0INT(calling, LIKELIHOOD_IGNORETREND); 
  // to do : betas_separate bereinigen && ignore_trend nicht alles
  // allocieren!

  // PMI0(calling);
  // printf("%d %d \n", LIKELIHOOD_BETASSEPARATE, PARAM0INT(calling, LIKELIHOOD_BETASSEPARATE));
  
  L->betas_separate = 
    CALLINGNR == LIKELIHOOD_CALL && PARAM0INT(calling, LIKELIHOOD_BETASSEPARATE);
  if (L->betas_separate && sets > 1) SERR("separate estimation of the betas only for one data set possible.");

  L->nas_boxcox =  L->nas_boxcox_mu = 0;
  for (i=0; i<two_vdim; i++) {
    L->nas_boxcox += ISNA(P(GAUSS_BOXCOX)[i]);
    if (i % 2 == 1) L->nas_boxcox_mu += ISNA(P(GAUSS_BOXCOX)[i]);
  }

  // information about free variables: 
  // first: get effect and NAs   
  
  assert(cov->calling != NULL && cov->calling->calling == NULL);
  ONCE_NEW_COV_STORAGE(cov->calling, mle);
  err = internalSetAndGet(cov, global->fit.lengthshortname,
			  LIKELI_NA_INTEGER, LIKELI_EXCLUDE_TREND,
			  LocxdimOZ(cov), global_var, facts,
			  original_model);
  if (err != NOERROR) goto ErrorHandling;

  // second: existence and "size" of known trend and fixed trend

  //  printf("\n\n\n\n .>>>>>>>> get Beta:\n");
  GetBeta(cov, L, &ne);
  // printf("\n\n\n\n .>>>>>>>> get Beta fodnex\n");

  if (L->fixedtrends + 1 < MAX_LIN_COMP)
    L->cum_n_betas[L->fixedtrends + 1] = L->cum_n_betas[L->fixedtrends]; 
  betatot = L->cum_n_betas[L->fixedtrends];

  //  printf("betatot = %d %d %d\n", betatot,L->dettrends,  L->fixedtrends);
  
  if (betatot > 0) {
    ne = 0;
    double **pt = L->where_fixed = 
      (double **) MALLOC( L->cum_n_betas[L->fixedtrends] * sizeof(double**) );
    GetBeta(cov, L, &ne, &pt);  // pt wird hochgezaehlt
    //    printf("%d %d\n", betatot, pt - L->where_fixed );
    assert(pt == L->where_fixed + betatot);
  }
 
  if (n < L->maxbeta) n = L->maxbeta;

  //printf("maxbeta %d %d %d\n", L->maxbeta, vdim, n); 

  // alloc_cov is need for both CovarianceMatrix and FctnInternXal
  // APMI0(cov);
  if ((err = alloc_fctn(cov, OWNTOTALXDIM, n)) !=NOERROR) goto ErrorHandling;
  
  //  printf("betas = %d/%d %d;  %d %d %d\n", L->cum_n_betas[ne], betatot, facts->neffect, facts->effect[0], facts->effect[1], facts->effect[2]);
  for(k=0; k<facts->neffect; k++) {
    //    printf("k=%d %d \n", k, facts->effect[k]);
    if (facts->effect[k] == DetTrendEffect) {
      L->nas_det[dummy_det++] = facts->nas[k];
    } else if (facts->effect[k] == FixedEffect) {
      L->nas_fixed[dummy_fixed++] = facts->nas[k];
    } 
  }

  assert(dummy_det == L->dettrends && dummy_fixed == L->fixedtrends && 
	 0 == L->random);

  for (int set = 0; set < sets; set++){
    cov->base->set = set;
    int totptsvdim = LoctotalpointsY(cov) * vdim;
    if (totptsvdim > L->max_total_data) 
      L->max_total_data = totptsvdim;
  }
  int Sets; Sets = sets; // 25.1.21 to fool the compiler
  L->X = (double **) CALLOC(Sets, sizeof(double *));
  L->YhatWithoutNA = (double **) CALLOC(Sets, sizeof(double *));
  for (i=0; i<sets; i++) L->X[i] = L->YhatWithoutNA[i] = NULL;
  L->sumY = (double *) MALLOC(L->max_total_data * sizeof(double)); 
  L->sets = sets;
  
  for (int set = 0; set < sets; set++){
    cov->base->set = set;
    double 
      *val = L->sumY;	  // dummy  
    Long 
      totptsvdim = LoctotalpointsY(cov) * vdim;
    
    double *Yhat = L->YhatWithoutNA[set] =
      (double*) CALLOC(totptsvdim, sizeof(double));

    //PMI0(cov);
    
    for (i=0; i<L->dettrends; i++) {
      if (L->nas_det[i] == 0) {
  	FctnIntern(cov, L->det_effect[i], L->det_effect[i], false, val);   
	for (k = 0; k < totptsvdim; k++) Yhat[k] += val[k];
      } else L->dettrend_has_nas = true; 
    }

    //    for (k = 0; k < totptsvdim; k++) printf("%10g ", Yhat[k]); printf("\n"); BUG;
    
    if (L->random) 
      GERR("Currently the linear part of random effects cannot be determined");
  }
 
  if (CALLINGNR == LINEARPART_CALL) {
    if (L->fixedtrends) {
      for (int set = 0; set < sets; set++){
	cov->base->set = set;
	Long tot = (Long) betatot * LoctotalpointsY(cov) * vdim;
	L->X[set] = (double*) CALLOC(tot, sizeof(double)); 
      }  
    }
  } else {
    listcpy(&(L->datasets),  PARAMLIST(calling, LIKELIHOOD_DATA), false);
    L->data_nas = (int *) MALLOC(sets * sizeof(int));
    listoftype *datasets = L->datasets;
       
    for (int set = 0; set < sets; set++){
      cov->base->set = set;
     double 
	*data = SET_OUT_OF(datasets);	    
      int
	ncol = NCOL_OUT_OF(datasets),
	repet = ncol / vdim,
	nrow = NROW_OUT_OF(datasets),
	data_nas = 0,
	ndata = nrow * ncol;

      //PMI0(cov->calling);   printf("nrow = %d %d\n", nrow, ncol);
 
      if (nrow != LoctotalpointsY(cov)) BUG;
      if (repet > maxrepet) maxrepet = repet;
      if (nrow  > maxtotpts) maxtotpts = nrow;
      if (ndata > maxndata) maxndata = ndata;
          
      // any_data_nas = false;  
      for (k = 0; k<ndata; data_nas += ISNA(data[k++]));
      L->data_nas[set] = data_nas;
      L->data_has_nas |= data_nas > 0;

      int vdimMax = MIN(vdim, MAXBOXCOXVDIM),
	total = vdim * nrow;
      if (L->nas_boxcox_mu) {
	bool warn_boxcox = true;
	for (int r=0; r<repet; r++) {
	  int z = r * total;
	  for(i=0; i<vdimMax; i++) {
	    double negmu = - global->fit.BC_lambdaLB[2 * i + 1];
	    for(k=0; k<nrow; k++, z++) {
	      if (!ISNA(data[z]) && data[z] <= negmu) {
		if (warn_boxcox) {
		  warning("data values do not ensure positive values with given lower bounds for the second parameter in the Box-Cox transformation, hence the lower bound has been increased.");
		  warn_boxcox = false;
		}
		negmu = data[k] - 1e-10;
		global->fit.BC_lambdaLB[2 * i + 1] = -negmu;	  
	      }
	    } // k, nrow
	  } // i, vdim
	} // r, repet
      } else if (R_FINITE(P(GAUSS_BOXCOX)[0])) {
	bool warn_boxcox = false;
	for (int r=0; r<repet; r++) {
	  int z = r * total;
	  for(i=0; i<vdimMax; i++) {
	    double negmu = - global->fit.BC_lambdaLB[2 * i + 1];
	    for(k=0; k<nrow; k++, z++) {
	      if ((warn_boxcox = !ISNA(data[z]) && data[z] <= negmu)) break;
	    }
	  }
	}
	if (warn_boxcox) warning("data values do not ensure positive values in Box-Cox transformation, hence likely errors will occur.");
      }
    }
 
    // provide necessary space and account for known information
    max_total_data_Sq = (Long) L->max_total_data * L->max_total_data;
    Long totdata_bytes = (Long) max_total_data_Sq * sizeof(double);

    //    PMI(cov);
    //    printf("totdata_bytes = %ld %ld %ld\n", totdata_bytes, L->max_total_data, L->max_total_data);
    
    L->C = (double *) MALLOC(totdata_bytes);
    //    printf("L->C %d\n", totdata_bytes);
    if (L->data_has_nas) {
      L->Cwork = (double *) MALLOC(totdata_bytes);
      Long tot = (Long) L->max_total_data + betatot * LoctotalpointsY(cov);
      L->Xwork = (double*) MALLOC(tot * sizeof(double));
      //printf("XWORK tot = %d\n", tot);
    }
    if (L->fixedtrends) {
      L->XitXi =  (double*) MALLOC(sizeof(double) * betatot * 
				   (betatot > maxrepet ? betatot : maxrepet)); 
      L->XtX = (double*) CALLOC(betatot * betatot, sizeof(double)); 
      int all_betatot = betatot;
      if (L->betas_separate) all_betatot *= maxrepet;
      L->XCY = (double *) MALLOC(all_betatot * sizeof(double)); 
      L->betavec = (double *) MALLOC(all_betatot * sizeof(double)); 
    }

    L->CinvXY = (double*) MALLOC(sizeof(double) *
				 ((Long)betatot* maxtotpts* vdim + maxndata)); 
    
    for (int set = 0; set<sets; set++){
      cov->base->set = set;
      double 
	*data = SET_OUT_OF(datasets),
	*YhatWithoutNA = L->YhatWithoutNA[set];	    
      int  m, p,
	ncol = NCOL_OUT_OF(datasets),
	repet = ncol / vdim,
	nrow = NROW_OUT_OF(datasets),
	ndata = nrow * ncol,
	totptsvdim = nrow * vdim;
      
      L->X[set] = 
	(double*) CALLOC((Long)(betatot + repet) * vdim * LoctotalpointsY(cov), 
                           sizeof(double)); // auch die Y-Werte; vdim included!
      
      
      if (!L->dettrend_has_nas && L->nas_boxcox == 0) {
	double *Xdata = L->X[set] + betatot * totptsvdim;
	//	printf("%d r=%d c=%d v=%d t=%d %d %d beta=%d repet=%d pts=%d %s\n",
	//	       ndata, nrow, ncol, vdim, totptsvdim,
	//       (betatot + repet * vdim) * LoctotalpointsY(cov),
		 //       betatot * totptsvdim, betatot, repet,LoctotalpointsY(cov),
	//	       NICK(cov));
	MEMCOPY(Xdata, data, ndata * sizeof(double));

	//  printf("data=%10g %d %d %d\n", data[0], ndata, betatot, totptsvdim);	
	//  BUG;
 
	for (p = m = 0; m < repet; m++) {
	  for (k = 0; k < totptsvdim; k++) {
	    Xdata[p++] -= YhatWithoutNA[k]; // XX lesefehler
	  }
	}
      
	if (L->nas_boxcox == 0)
	  boxcox_trafo(P(GAUSS_BOXCOX), vdim, Xdata, nrow, repet);    
	
      }

      for (i=0; i<L->random; i++) {
	BUG; // not programmed yet
	if (L->nas_random[i] == 0) {
	  BUG;
	} else L->random_has_nas = true;
      }
    } // global->set
  } // not rftrend

  for (int set = 0; set < sets; set++){
    cov->base->set = set;
    for (i=0; i<L->fixedtrends; i++) {
      if (L->nas_fixed[i] == 0) {
	FctnIntern(cov, L->fixed_effect[i], L->fixed_effect[i], false,
		   L->X[set] + L->cum_n_betas[i] * LoctotalpointsY(cov));// "Lueckentext"!!
	//	printf("hier %d %d incr=%d %10g\n", L->fixedtrends, L->nas_fixed[i],
	//	       L->cum_n_betas[i] * LoctotalpointsY(cov),
	//	       L->X[set] [L->cum_n_betas[i] * LoctotalpointsY(cov)]);
	//	APMI(cov);
	
      }
      else L->fixedtrend_has_nas = true;
    }

    //    int kk = 0;
    //    for (int h=0; h<betatot; h++) { printf("\n");
    //      for (int hh=0; hh<LoctotalpointsY(cov); hh++) {
    //	printf("%d ", (int)L->X[set][kk++]);
    //      }
    //    } BUG;

  } // global->set

  EXT_NEW_STORAGE(solve);
  

  //  assert(cov->Slikelihood != NULL);  printf("%d\n", addressbits(cov->Slikelihood));  APMI(cov);

 ErrorHandling:
  cov->base->set = 0;
  RETURN_ERR(err);
} // struct logli
  
    
