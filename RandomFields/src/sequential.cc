/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by sequential method

 Copyright (C) 2001 -- 2017 Martin Schlather, 

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
#include <stdio.h>  
//#include <stdlib.h>
#include <R_ext/Lapack.h>
#include "questions.h"
#include "Processes.h"
#include "Coordinate_systems.h"


//#define  debug_sequ 1

bool debugging = true;



#define SEQU_BACK (COMMON_GAUSS + 1)
#define SEQU_INIT (COMMON_GAUSS + 2)

static inline double scalar(double *A, double *B, Long C) {
 double  sum=0.0;			    
 for (Long i=0; i<C; i++) sum += A[i] * B[i];
 return sum;
}
//#define SCALAR_PROD(A,B,C) scalar(A,B,C)
#define SCALAR_PROD(A,B,C) Ext_scalarX(A,B,C, SCALAR_AVX)
 


int check_sequential(model *cov) {
  globalparam *global = &(cov->base->global);
#define nsel 4
  model *next=cov->sub[0];
  int err,
    dim = ANYDIM; // taken[MAX DIM],
  sequ_param *gp  = &(global->sequ);

  if (!Locgrid(cov) && !LocTime(cov)) 
    SERR1("'%.50s' only possible if at least one direction is a grid", NICK(cov));

  kdefault(cov, SEQU_BACK, gp->back);
  kdefault(cov, SEQU_INIT, gp->initial);
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  RESERVE_BOXCOX;

  
  if ((err = CHECK(next, dim, dim, PosDefType, XONLY, 
		   SymmetricOf(OWNISO(0)), // todo : eigentlich sollte coordinate of reichen; aber unten alg durchschauen.
		   SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
    RETURN_ERR(err);
  }
  
  if (next->pref[Sequential] == PREF_NONE) RETURN_ERR(ERRORPREFNONE);
  setbackward(cov, next);

  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov)) != NOERROR) RETURN_ERR(err);


  RETURN_NOERROR;
}


void range_sequential(model VARIABLE_IS_NOT_USED *cov, int k, int i, int j,
		       simple_range_type *range) {
  switch(k) {
    GAUSS_BOXCOX_RANGE;
  case SEQU_BACK:
    range->min = 0;
    range->max = RF_INF;
    range->pmin = 0.1;
    range->pmax = 10;
    range->openmin = true;
    range->openmax = true; 
    break;
  case SEQU_INIT:
    range->min = RF_NEGINF;
    range->max = RF_INF;
    range->pmin = - 10;
    range->pmax = 10;
    range->openmin = false;
    range->openmax = true;
    break;
  default : BUG;  
  }
}

// (U1, U2) ~ N(0, S) =>
// U1 | U2 ~ N(S_{12} S_{22}^{-1} u2, S_11 - S12 S_22^{-1} S_21)

// start mit S_22; dann nur zeilenweise + Einschwingen

int init_sequential(model *cov, gen_storage VARIABLE_IS_NOT_USED *s){
  globalparam *global = &(cov->base->global);
  model *next = cov->sub[0];
  location_type *loc = Loc(cov);
  if (loc->distances) RETURN_ERR(ERRORFAILED);

  int withoutlast,  endfor,
    err = NOERROR,
    dim = ANYDIM,
    spatialdim = dim - 1,
    vdim = next->vdim[0],
    max = global->direct.maxvariables,
    back= P0INT(SEQU_BACK), 
    initial= P0INT(SEQU_INIT);
  assert(dim == loc->timespacedim);
  assert(next->vdim[0] == next->vdim[1]);
  if (initial < 0) initial = back - initial;

  double 
    //*caniso = loc->caniso,
    *timecomp = loc->grid ? loc->xgr[spatialdim] : loc->T,
    *xx = NULL,
    *G = NULL, 
    *U11 = NULL, 
    *COV21 = NULL,
    *MuT = NULL, 
    *U22 = NULL, 
    *Inv22 = NULL;
  double
    *res0 = NULL;
  sequ_storage* S = NULL;
  Long  
    timelength = (int) (loc->grid ? loc->xgr[spatialdim] :loc->T)[XLENGTH],
    spatialpnts = loc->totalpoints / timelength,
    totpnts = back * spatialpnts, 
    totpntsSQ =  totpnts * totpnts,
    spatialpntsSQ = spatialpnts * spatialpnts,
    spatialpntsSQback = totpnts * spatialpnts;

  bool 
    storing = cov->base->stored_init,
    Time = loc->Time;
   DEFAULT_INFO(info);

  if (hasEvaluationFrame(cov)) {
    RETURN_NOERROR;
  }
  
  cov->method = Sequential;
     
  if (!loc->grid && !Time) 
    GERR("last component must be truely given by a non-trivial grid");

  if (DefList[NEXTNR].implemented[Sequential] != IMPLEMENTED) {
    err=ERRORNOTDEFINED; 
    goto ErrorHandling;
  }
  
  if (VDIM0 > 1) {
      err=ERRORNOMULTIVARIATE; 
      goto ErrorHandling;   
  }
  
  if (totpnts > max) // (totpnts * vdim > max)
    GERR6("'%.50s' valid only if the number of locations is less than '%.50s' (=%d) . Got %d * %ld = %ld.", NICK(cov), direct[DIRECT_MAXVAR_PARAM],
	  max, back, spatialpnts, totpnts);
   
  if (timelength <= back) {
    GERR2("the grid in the last direction is too small; use method '%.50s' instead of '%.50s'",
	  DefList[DIRECT].nick, DefList[SEQUENTIAL].nick);
  } 
  if (back < 1) back = max / spatialpnts;

  //PMI0(cov);
  //  TREE0(cov);
  //  printf("spectal %d %d %d %d vdim=%d back=%d spa=%d tot=%d\n", totpntsSQ, spatialpntsSQback, totpnts, vdim * (totpnts + spatialpnts * initial), vdim, back, spatialpnts, totpnts);

  if ((U22 = (double *) MALLOC(sizeof(double) * totpntsSQ))==NULL ||
      (Inv22 = (double *) MALLOC(sizeof(double) * totpntsSQ))==NULL ||
      (COV21 = (double *) MALLOC(sizeof(double) * spatialpntsSQback))==NULL ||
      (U11 = (double *) MALLOC(sizeof(double) * spatialpntsSQ))==NULL ||
      (MuT = (double *) MALLOC(sizeof(double) * spatialpntsSQback))==NULL ||
      (G = (double *) MALLOC(sizeof(double) * totpnts))==NULL ||
      (res0 = (double *) MALLOC(sizeof(double) * vdim *
				  (totpnts + spatialpnts * initial))) ==NULL) {
    err=ERRORMEMORYALLOCATION;  
    goto ErrorHandling;
  }
  NEW_STORAGE(sequ);
  S = cov->Ssequ;

 
  /* ************************* */
  /* create matrix of explicit */
  /*       x-coordinates       */
  /* ************************* */

 
  if (loc->grid) loc->xgr[spatialdim][XLENGTH] = back; 
  else loc->T[XLENGTH] = back;
  TransformLoc(cov, &xx, DOLLAR_IMPOSSIBLE);
  loc = Loc(cov);
  assert(loc->caniso == NULL);
  if (loc->grid) loc->xgr[spatialdim][XLENGTH]=timelength; 
  else loc->T[XLENGTH] = timelength;
 
  /* ********************* */
  /* matrix creation part  */
  /* ********************* */
  int  row, sub_err;
  double *y;
  y = (double*) MALLOC(dim * sizeof(double));

  if (PL >= PL_SUBDETAILS) { LPRINT("covariance matrix...\n"); }
  for (Long k0=0, k=0, icol = 0; icol<totpnts; icol++, k0+=dim) {
    k += icol;
    for (Long segment = icol * (totpnts + 1), irow=icol, k2=k0; 
	 irow < totpnts; irow++) {
      Long k1 = k0;
      for (int d=0; d<dim; d++) {
	y[d] = xx[k1++] - xx[k2++];
      }
      COV(y, info, next, U22 + segment); 
      U22[k++] = U22[segment];
      segment += totpnts;
    }
  }

  /*
  for (k=i=0; i<totpnts; i++) {
    for (d=0; d<totpnts; d++) {
      printf("%10g ", U22[k++]); //
    }
//    printf("\n");
  }
  */
       
  /* Orndung von U22: 
      back ....1
   back
   ...
   1

  Ordnung von COV21
  back+1
  back
  ...
  2
  */

  // *** S11 und S21
  withoutlast = totpnts - spatialpnts;
  for (Long jj=0, k=0, i=withoutlast * totpnts; i<totpntsSQ; ) {
    jj += spatialpnts;
    endfor = jj + withoutlast;
    for (; jj<endfor; ) {
	COV21[jj++] = U22[i++];
    }
    endfor = k + spatialpnts;
    for (; k<endfor; ) U11[k++] = U22[i++];
  }
  

  
  // *** S21 rest
  y[spatialdim] = timecomp[XSTEP] * back;
  for (Long k=0, k0 = 0, icol = 0; icol<spatialpnts; icol++, k0+=dim) {//t_{n+1}
    for (Long k2=0, irow=0; irow<spatialpnts; irow++) { // t_1
	Long k1 = k0;
	for (int d=0; d<spatialdim; d++) {
	    y[d] = xx[k1++] - xx[k2++];
	}
	COV(y, info, next, COV21 + k);	
	k++; k2++;
    }
    k += withoutlast;
  }

  FREE(y);

//  for (i=0; i<withoutlast * spatialpnts; i++) {
//    print("%3.4f ", COV21[i]);
//  }
//  assert(false);


 #ifdef debug_sequ
  LPRINT("U11 spatialpts=%d\n", spatialpnts);
  for (int i=0; i<spatialpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
    }
    LPRINT("\n");
  }

  LPRINT("S22 %d %d\n", totpnts,spatialpnts );
  for (int i=0; i<totpnts; i++) {
    LPRINT("%2d: ", i);
    for (int j=0; j<totpnts; j++) {
      LPRINT("%3.3f ", U22[i + j * totpnts]);
   }
    LPRINT("\n");
  }
#endif



  /* ********************* */
  /* matrix compositions  */
  /* ********************* */

  // *** sqrt S 22 
  if (PL>=PL_STRUCTURE) { LPRINT("Cholesky decomposition of Cov22...\n"); }
  row=totpnts; 
  if ((sub_err = Ext_chol(U22, row)) != NOERROR) {
    if (PL>=PL_ERRORS) {
      LPRINT("Error code in cholesky = %d\n", sub_err);
    }
    err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  } 

   
  // *** inverse of S 22 
  if (PL>=PL_STRUCTURE) { LPRINT("inverse of Cov22...\n"); }
  MEMCOPY(Inv22, U22, sizeof(double) * totpntsSQ);
  Ext_chol2inv(Inv22, row);
  for (Long k=0, i=0; i<totpnts; i++) {
    k += i;
    for (Long segment=i* totpnts+i; segment<totpntsSQ; segment += totpnts) {
	Inv22[k++] = Inv22[segment]; 
//	LPRINT("sg %d %10g %d \n", k, Inv22[segment] , segment);
    }
  }



#ifdef debug_sequ
 LPRINT("Cholesky S22 %d %d\n", totpnts,spatialpnts );
  for (int i=0; i<totpnts; i++) {
    LPRINT("%2d: ", i);
    for (int j=0; j<totpnts; j++) {
      LPRINT("%3.3f ", U22[i + j * totpnts]);
   }
    LPRINT("\n");
  }
  
   LPRINT("Inverse S22 %d %d\n", totpnts,spatialpnts );
   for (int i=0; i<totpnts; i++) {
    LPRINT("%2d: ", i);
    for (int j=0; j<totpnts; j++) {
      LPRINT("%3.3f ", Inv22[i + j * totpnts]);
   }
    LPRINT("\n");
  }
  
  
  LPRINT("C21\n");
  for (int i=0; i<totpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.3f ", COV21[i + j * totpnts]);
   }
    LPRINT("\n");
  }

#endif


  // *** Mu = C21^T Inv22
  // b)ack s)spatial
  // MuT : sb x s
  // Mu : s x sb
  // C21 : sb x s
  // Inv22 : (sb)^2
  if (PL>=PL_STRUCTURE) { LPRINT("calculating matrix for the mean...\n"); }

  // ##define remove_bottom_zeros 1
#ifdef remove_bottom_zeros  
  double *Maxi;
  Maxi = (double*) CALLOC(back, sizeof(double));
#endif  
  for (Long l=0, z=0; z<spatialpnts; z++) { // C21:s
    int i = z * totpnts;
    for (int k=0, j=0; k<back; k++) { // Inv sb
      double m=RF_NEGINF;
      for (int n=0; n<spatialpnts; n++, j+=totpnts, l++) { // Inv sb
	double dummy = 0.0;
	// for (k=0; k<totpnts; k++) dummy += COV21[i + k] * Inv22[j + k];
	//	MuT[l++] = dummy;
	MuT[l] = SCALAR_PROD(COV21 + i, Inv22 + j, totpnts);
#ifdef remove_bottom_zeros  	
	double f = FABS(MuT[l]);
	m = MAX(m, f);
#endif  
      }
#ifdef remove_bottom_zeros  
      Maxi[k] = MAX(Maxi[k], m);
 #endif  
   }
  }


  S->delta_back_MuT = 0;
  
#ifdef remove_bottom_zeros  
  double threshold;
  threshold = 0.0;
  for (int i=0; i<back; i++) threshold = MAX(threshold, Maxi[i]);
  threshold *= 1e-13;
  // nett gedacht mit den Nullen wegstreichen. Bringt aber quasi nichts,
  // nicht mal bei exponential
  for (int start=0; start<back; start++) {
     if (Maxi[start] > threshold) {
      S->delta_back_MuT = start;
      // printf("gain = %d\n",  S->delta_back_MuT);
      break;
    }
  }
#endif  
 


#ifdef debug_sequ
  LPRINT("MuT\n");
  for (int i=0; i<totpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2e", MuT[i + j * totpnts]);
   }
    LPRINT("\n");
  }
#endif

  //  assert(false);
 
  // *** C11
  if (PL>=PL_STRUCTURE) { LPRINT("calculating cov matrix...\n"); }
  for (Long l=0, i=0; i<spatialpntsSQback; i+=totpnts) {
    for (int j=0; j<spatialpntsSQback; j+=totpnts) {
      // TO DO: SYMMETRIE NUTZEN
	double dummy = 0.0;
	for (int kk=0; kk<totpnts; kk++) dummy += MuT[i + kk] * COV21[j + kk];
	U11[l++] -= SCALAR_PROD(MuT + i, COV21 + j, totpnts);
    }
  }

/*
  LPRINT("U11\n");
  for (i=0; i<spatialpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
   }
   LPRINT("\n");
  }
*/


  // *** sqrt C11
  if (PL>=PL_STRUCTURE) { LPRINT("Cholesky decomposition of Cov11...\n"); }
  row=spatialpnts;
  
  sub_err = Ext_chol(U11, row);

#ifdef debug_sequ
  LPRINT("trafo U11 spatialpts=%d\n", spatialpnts);
  for (int i=0; i<spatialpnts; i++) {
    for (int j=0; j<spatialpnts; j++) {
      LPRINT("%3.2f ", U11[i + j * spatialpnts]);
    }
    LPRINT("\n");
  }
#endif
  
    if (sub_err!=NOERROR) {
    if (PL>=PL_ERRORS) { LPRINT("U11: Error code in chol = %d\n", sub_err); }
    err=ERRORDECOMPOSITION;
    goto ErrorHandling;
  }

  err = ReturnOwnField(cov);

 ErrorHandling: // and NOERROR...
  FREE(xx);
  if (S!=NULL) {
      S->totpnts = totpnts;
      S->spatialpnts = spatialpnts;
      S->back = back;
      S->initial = initial;
      S->ntime = (int) timecomp[XLENGTH];
  }
  if (!storing && err!=NOERROR) {
    //  printf("storing = %d\n", storing); BUG;
    FREE(MuT);
    FREE(U11);
    FREE(U22);
    FREE(G); 
    FREE(res0); 
  } else {
    if (S != NULL) {
      // printf("%f %f\n", U11[0], U11[1]);
      S->U22=U22;
      S->U11=U11;
      S->MuT=MuT;
      S->G=G;
      S->res0=res0;
      //assert(false);
     }
  }
  if (COV21!=NULL) {
      if (S != NULL && debugging)  S->Cov21 = COV21;
      else UNCONDFREE(COV21);
  }
  if (Inv22!=NULL) {
      if (S != NULL && debugging)  S->Inv22 = Inv22;
      else UNCONDFREE(Inv22);
  }

  cov->simu.active = err == NOERROR;

  RETURN_ERR(err);
}


void sequentialpart(double *res, Long totpnts, int spatialpnts,
		    int truelybackpnts, // == totpnts by default
		    int ntime,
		    double *U11, double *MuT, double *G) {
  MuT += totpnts - truelybackpnts;
  double *rp = res + totpnts;
  for (int n=0; n<ntime; n++, rp += spatialpnts) {    
    for (int i=0; i<spatialpnts; i++) G[i] =GAUSS_RANDOM(1.0);
#ifdef DO_PARALLEL
#pragma omp parallel for num_threads(CORES) if (spatialpnts > 5 * CORES)
#endif
    for (int i=0; i<spatialpnts; i++) {
      rp[i] = SCALAR_PROD(G, U11 + i * spatialpnts, i+1) 
	+ SCALAR_PROD(MuT + i * totpnts, rp - truelybackpnts, truelybackpnts);
    }
  }
}


void do_sequential(model *cov, gen_storage VARIABLE_IS_NOT_USED *s) 
{  
  model *next = cov->sub[0];
  sequ_storage
    *S = cov->Ssequ;
  assert(S != NULL); 
  
  int vdim = next->vdim[0];
  Long totpnts = S->totpnts;
  double *G,*U22, *U11, *MuT;
  double *res0,
    *res = cov->rf; 
  SAVE_GAUSS_TRAFO;
 
  assert(res != NULL);
  assert(S != NULL); 

 
  //  printf("totpnts %ld %d %d\n", totpnts, S->initial, S->spatialpnts);
  // assert(false);
 
  U22 = S->U22;  // S22^{1/2}
  U11 = S->U11;  // S11^{1/2}
  MuT = S->MuT;
  res0 = S->res0;
  G = S->G;// only the memory space is of interest (stored to avoid 
  //          allocation errors here)

  // Start
  for (int i=0; i<totpnts; i++) G[i] = GAUSS_RANDOM(1.0);
  for (int k=0, i=0; i<totpnts; i++, k+=totpnts){
    double dummy, *Uk;
    Uk = &U22[k]; 
    dummy =0.0;
    for (int j=0; j<=i; j++){
      // printf("k=%d %d %f %f\n", k, j, G[j], Uk[j]);
      dummy += G[j] * Uk[j];
    }
    res0[i] = (double) dummy;
    // printf("%d %f \n",i, res0[i]);
  }

  int truely_back_points = totpnts - S->spatialpnts * S->delta_back_MuT;
  sequentialpart(res0, totpnts, S->spatialpnts, truely_back_points, 
		 S->initial, U11, MuT, G);
  res0 += S->initial * S->spatialpnts;
  MEMCOPY(res, res0, sizeof(double) * totpnts * vdim);

  sequentialpart(res, totpnts, S->spatialpnts, truely_back_points,
		 S->ntime - S->back, U11, MuT, G);
 

  BOXCOX_INVERSE;
}
