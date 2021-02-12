
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Simulation of a random field by circulant embedding
 (see Wood and Chan, or Dietrich and Newsam for the theory)
 and variants by Stein and by Gneiting

 Copyright (C) 2001 -- 2003 Martin Schlather
 Copyright (C) 2004 -- 2005 Yindeng Jiang & Martin Schlather
 Copyright (C) 2006 -- 2011 Martin Schlather
 Copyright (C) 2011 -- 2015 Peter Menck & Martin Schlather

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

#include <R_ext/Linpack.h>
#include <R_ext/Utils.h>     
#include <R_ext/Lapack.h> // MULT 

//#include <stdlib.h>
#include "def.h"

#include "questions.h"
#include "operator.h"
#include "Processes.h"


#ifdef CE_PARALLEL
#include <omp.h>
#endif






/*********************************************************************/
/*           CIRCULANT EMBEDDING METHOD (1994) ALGORITHM             */
/*  (it will always be refered to the paper of Wood & Chan 1994)     */
/*********************************************************************/


int fastfourierInit(int *m, int dim, FFT_storage *FFT)
  /* this function is taken from the fft function by Robert Gentleman 
     and Ross Ihaka, in R */
{
  int Err, maxp, maxmaxf,maxmaxp,
    nseg = maxmaxf =  maxmaxp = 1;
  /* do whole loop just for err checking and maxmax[fp] .. */
  
 
  for (int i = 0; i<dim; i++) {
    if (m[i] > 1) {
      if (fft_factor(m[i], FFT->maxf + i, &maxp, FFT->kt + i, FFT->m_fac + i,
		     FFT->NFAC[i]))
	FAILED("fft factorization failed");
      
      if (FFT->maxf[i] > maxmaxf) maxmaxf = FFT->maxf[i];
      if (maxp > maxmaxp) maxmaxp = maxp;
      nseg *= m[i];
    }
  }
  
  FREE(FFT->work);
  FREE(FFT->iwork);
  //
  if ((FFT->work = (double*) MALLOC(4 * maxmaxf * sizeof(double)))==NULL || 
	(FFT->iwork = (int*) MALLOC(maxmaxp  * sizeof(int)))==NULL) { 
    Err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  FFT->nseg = nseg;

  //  nseg = leng th(z); see loop above
  return NOERROR;

ErrorHandling:
  return Err;
}


int fastfourier(double *data, int *m, int dim, bool first, bool inverse,
		FFT_storage *FFT) {
  int Err;
  if (first && (Err = fastfourierInit(m, dim, FFT)) != NOERROR) return Err;
  return fastfourier(data, m, dim, inverse, FFT);
}


int fastfourier(double *data, int *m, int dim,  bool inverse,
		FFT_storage *FFT) {
 /* this function is taken from the fft function by Robert Gentleman 
     and Ross Ihaka, in R */
  int nfac[21];
  Long
    inv = (inverse) ? 2 : -2,
    n = 1,
    nspn = 1,
    nseg = FFT->nseg;
  for (Long i = 0; i < dim; i++) {
    if (m[i] > 1) {
      nspn *= n; 
      n = m[i]; 
      nseg /= n;
      MEMCOPY(nfac, FFT->NFAC[i], sizeof(int) * 21);
      //     printf("%d: %d %d %d; fnax=%d %d %d %d %d %d %d %d \n", L, FFT->maxf[i], FFT->kt[i], FFT->m_fac[i], nfac[0], nfac[1], nfac[2], nfac[3], nfac[4], nfac[5], nfac[6], nfac[7]);
      if (!fft_work(data, data + 1, nseg, n, nspn, inv, FFT->work, FFT->iwork,
		    FFT->maxf[i], FFT->kt[i], FFT->m_fac[i], nfac))
	return XERRORFOURIER;
    }
  }
  return NOERROR;  
}


#define LOCPROC_DIAM  (CE_LAST + 1)// parameter p
#define LOCPROC_A  (CE_LAST + 2) // always last of LOCPROC_*
#define LOCPROC_R LOCPROC_A // always last of LOCPROC_*


int init_circ_embed(model *cov, gen_storage VARIABLE_IS_NOT_USED  *S){
   globalparam *global = &(cov->base->global);
 int err=NOERROR;
  if (!Locgrid(cov)) SERR("circ embed requires a grid");

  double  steps[MAXCEDIM], 
      **c = NULL, // For multivariate *c->**c
      **Lambda = NULL;
  Long  hilfsm[MAXCEDIM],
    mtot=0;
  model *next = cov->sub[0];
  location_type *loc = Loc(cov);
  double
    maxGB = P0(CE_MAXGB),
    tolRe = P0(CE_TOLRE),
    tolIm = P0(CE_TOLIM),
    *caniso = loc->caniso,
    *mmin = P(CE_MMIN);
  int dim = OWNLOGDIM(0),
    dimM1 = dim - 1,
    vdim   = VDIM0,
    cani_ncol = loc->cani_ncol,
    cani_nrow = loc->cani_nrow,
    lenmmin = cov->nrow[CE_MMIN],
    vdimSq = vdim * vdim, // PM 12/12/2008
    maxmem = P0INT(CE_MAXMEM),
    trials = global->internal.examples_reduced ? 1 : P0INT(CE_TRIALS),
    strategy = P0INT(CE_STRATEGY),
    worksize = vdim * 33;
  bool 
    isaniso = loc->caniso != NULL,
    force = (bool) P0INT(CE_FORCE),
    useprimes = (bool) P0INT(CE_USEPRIMES),
    dependent = (bool) P0INT(CE_DEPENDENT);
  double *tmp = (double *) MALLOC(vdim * vdim * sizeof(double)),	
    *tmpLambda = (double *) MALLOC(vdim * sizeof(double)), 
    *rwork = (double *) MALLOC( (vdim * 3 -2) * sizeof(double));
  complex
    *R = (complex*) MALLOC(vdim * vdim * sizeof(complex)),
    *work = (complex *) MALLOC(worksize * sizeof(complex));
  DEFAULT_INFO(info);

 
assert(VDIM0 == VDIM1);
  if (loc->distances) RETURN_ERR(ERRORFAILED);
  if (dim > MAXCEDIM) RETURN_ERR(ERRORMAXDIMMETH);
 
  cov->method = CircEmbed;

  NEW_STORAGE(ce);
  getStorage(s ,   ce); 
#ifdef CE_PARALLEL
  //  printf("%d %d\n", vdimSq, sizeof(FFT_storage));
  s->FFT = (FFT_storage *) CALLOC(vdimSq,  sizeof(FFT_storage)); 
#endif  
  int *mm = s->m,
    *halfm= s->halfm,
    *cumm = s->cumm;
  // since 20 Aug 17:  in NEW_STORAGE
  //  s->trials=0;
  //  s->cur_call_odd = 0;  // first one is set to be odd at beginning of
  //                        // do_circ_embed
  //  s->positivedefinite = false;     
  /* Eq. (3.12) shows that only j\in I(m) [cf. (3.2)] is needed,
     so only the first two rows of (3.9) (without the taking the
     modulus of h in the first row)
     The following variable 'index' corresponds to h(l) in the following
way: index[l]=h[l]        if 0<=h[l]<=mm[l]/2
index[l]=h[l]-mm[l]   if mm[l]/2+1<=h[l]<=mm[l]-1     
Then h[l]=(index[l]+mm[l]) % mm[l] !!
   */

  for (int d=0; d<dim; d++) {
    s->nn[d]=(int) loc->xgr[d][XLENGTH]; 
    steps[d]=loc->xgr[d][XSTEP];
  }
  s->vdim = vdim;

  // PM 12/12/2008
  // in case of multiple calls: possible to use imaginary part of
  // field generated in preceding run. If s->dependent is set, fields
  // are read off in the following order:
  //  1. real field 1st square
  //  2. imag field 1st square
  //  3. real field 2nd square
  //  4. imag field 2nd square
  //  ...
  //  2k+1. real field (2k+1)th square
  //  2k+2. imag field (2k+1)th square
  //
  //  i.e. read off imag fields in even calls. Generate new ones in
  //  odd calls resp. read off real squares if s->dependent==1



  /* cumm[i+1]=\prod_{j=0}^i m[j] 
     cumm is used for fast transforming the matrix indices into an  
     index of the vector (the way the matrix is stored) corresponding 
     to the matrix                   */                               
  /* calculate the dimensions of the matrix C, eq. (2.2) in W&C */       

  //  ("CE missing strategy for matrix entension in case of anisotropic fields %d\n)
  //  for (i=0;i<dim;i++) { // size of matrix at the beginning  
  //	print("%d %10g\n", i, mmin[i]);
  //  }

  bool multivariate;
  multivariate = VDIM0 > 1;
  assert(VDIM0 == VDIM1);

  // multivariate = !true; // checken !!

  for (int i=0;i<dim;i++) { // size of matrix at the beginning  
    double hilfsm_d;
    if (FABS(mmin[i % lenmmin]) > 10000.0) {     
      GERR1("maximimal modulus of mmin is 10000. Got %10g", mmin[i % lenmmin]);
    }
    hilfsm_d = (double) s->nn[i];
    if (mmin[i % lenmmin]>0.0) {
      if (hilfsm_d > (1 + (int) CEIL(mmin[i % lenmmin])) / 2) { // plus 1 since 
	// mmin might be odd; so the next even number should be used
        GERR3("Minimum size in direction %d is %10g. Got %10g\n",
	      (int) i, hilfsm_d, CEIL(mmin[i % lenmmin]));
      }
      hilfsm_d = (double) ( (1 + (int) CEIL(mmin[i % lenmmin])) / 2);
    } else if (mmin[i % lenmmin] < 0.0) {
      assert(mmin[i % lenmmin] <= - 1.0);
      hilfsm_d = CEIL((double) hilfsm_d * - mmin[i % lenmmin]);
    }
    if (hilfsm_d >=  MAXINT) {err=ERRORMEMORYALLOCATION; goto ErrorHandling;}
    if (useprimes) {
      // Note! algorithm fails if hilfsm_d is not a multiple of 2
      //       2 may not be put into NiceFFTNumber since it does not
      //       guarantee that the result is even even if the input is even !
	 hilfsm_d = 2.0 * NiceFFTNumber((int) hilfsm_d);
    } else {
      double dummy;
      if (multivariate)
	dummy = POW(3.0, 1.0 + CEIL(LOG(hilfsm_d) / LOG3 - EPSILON1000));
      else
	dummy = POW(2.0, 1.0 + CEIL(LOG(hilfsm_d) * INVLOG2-EPSILON1000));
      hilfsm_d =ROUND(dummy);
    }
    if (hilfsm_d >=  MAXINT) {err=ERRORMEMORYALLOCATION; goto ErrorHandling;}
    hilfsm[i] = (int) hilfsm_d;
  }
   

  
  // print("trials %d strat=%d force=%d prim=%d dep=%d  max=%10g re=%10g im=%10g\n", 
  //	 trials, strategy, force, useprimes, dependent, 
  //	 maxmem, tol_re, tol_im);


 

  /* The algorithm below is as follows:
     while (!s->positivedefinite && (s->trials<trials)){
     (s->trials)++;
     calculate the covariance values "c" according to the given "m"
     fastfourier(c)
     if (!force || (s->trials<trials)) {
     check if positive definite
     if (!s->positivedefinite && (s->trials<trials)) {
     enlarge "m"
     }
     } else 
     print "forced" //
     }
   */

  if (PL>=PL_STRUCTURE) { LPRINT("calculating the Fourier transform\n"); }
  while (!s->positivedefinite && (s->trials < trials)) {
    // print("trials %d %d\n", trials, s->trials);

    (s->trials)++;
    cumm[0]=1; 
    double realmtot = 1.0;
    for (int i=0; i<dim; i++) {
      if (hilfsm[i] < MAXINT) {
	mm[i] = hilfsm[i];
      } else {
	GERR3("Each direction allows for at most %d points. Got %d in the %d-th direction", (int) MAXINT, (int) hilfsm[i], i);
      }
      halfm[i]=mm[i] / 2; 
      cumm[i+1]=cumm[i] * mm[i]; // only relevant up to i<dim-1 !!
      realmtot *= (double) mm[i];
    }
    s->mtot = mtot = cumm[dim];

    if (PL>=PL_BRANCHING) {
      LPRINT("Memory need for ");
      for (int i=0; i<dim; i++) { 
	if (i > 0) { PRINTF(" x "); }
 	PRINTF("%d", mm[i]); 
      }
      PRINTF(" %s is 2 x %ld complex numbers.\n", dim == 1 ? "locations" : "grid", mtot  * vdimSq);
    }

    if ( (maxmem != NA_INTEGER && realmtot * vdimSq > maxmem))
      GERR4("total real numbers needed=%.0f > '%.50s'=%d -- increase '%.50s'",
	    realmtot * vdimSq, CE[CE_MAXMEM - COMMON_GAUSS - 1], 
	    maxmem, CE[CE_MAXMEM - COMMON_GAUSS - 1])
    else if (R_FINITE(maxGB) && realmtot * vdimSq * 32e-9 > maxGB)
      GERR4("total memory needed=%4.4f GB > '%.50s'=%4.4f GB -- increase '%.50s'",
	    realmtot* vdimSq * 32e-9, CE[CE_MAXGB - COMMON_GAUSS - 1], 
	    maxGB, CE[CE_MAXGB - COMMON_GAUSS - 1]);

    if (c != NULL) {
      for (int l=0; l<vdimSq; l++) FREE(c[l]);
      UNCONDFREE(c);
    }
    // for the following, see the paper by Wood and Chan!
    // meaning of following variable c, see eq. (3.8)
    if ((c = (double **) MALLOC(vdimSq * sizeof(double *))) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    for (int l=0; l<vdimSq; l++) {
      //     printf("%d maxmem=%d %d %10g\n", 2 * mtot * sizeof(double), maxmem, mtot, realmtot);
      if( (c[l] = (double *) MALLOC(2L * mtot * sizeof(double))) == NULL) {
	err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
      //printf("delte");      
      for (int k=0; k<2 * mtot; k++) c[l][k] = RF_NA;
    }
    
    Long err_occurred = NOERROR;
 #ifdef CE_PARALLEL
    //#pragma omp parallel for num_threads(CORES) reduction (+:err_occurred) if (vdim > 1L)
#endif
    
    for (int i=0; i<mtot; i++) {
      //     if (i==0) printf("A %d\n", omp_get_num_threads());
     int index[MAXCEDIM],
       dummy = i << 1;
     double hx[MAXCEDIM];
      bool cur_crit = false;
      SPLIT(i, mm, dim, index);

      for (int k=0; k<dim; k++) {
	int idxk = index[k];
	hx[k] = -steps[k] *// with '-' consistent with x[d]-y[d] in nonstat2stat
	  (double) ((idxk <= halfm[k]) ? idxk : idxk - mm[k]);
	cur_crit |= (idxk==halfm[k]);
      }
 
      //      printf("ci %d %10e %10e \n", dim, hx[0], hx[1]);

     if (isaniso) {
	double z[MAXCEDIM];
	xA_noomp(hx, caniso, cani_nrow, cani_ncol, z);
	COV(z, info, next, tmp);
     } else COV(hx, info, next, tmp);

     // printf("-- %d %10e %10e \n", dim, hx[0], hx[1]);
     //     printf("\n");
 
      for(int l=0; l<vdimSq; l++) {
	c[l][dummy] = tmp[l];
	c[l][dummy+1] = 0.0;
	
	if (ISNAN(c[l][dummy])) {

#ifdef CE_PARALLEL 
	  err_occurred++; break;
#else	     
	  GERR1("covariance models returns NAs. %.200s", CONTACT);
#endif
	}
	//	  print("%d %d %10g\n", i, mtot,  c[l][dummy]);
      }

      if (cur_crit) {
	for (int k=0; k<dim; k++) {
	  if (index[k]==halfm[k]) hx[k] -= steps[k] * (double) mm[k];
	}
	COV(hx, info, next, tmp);
	//	printf("c%10g %10g %10g %10g\n", tmp[0], tmp[1], tmp[2], tmp[3]);
	for (int l=0; l<vdimSq; l++) {
	  c[l][dummy] = 0.5 *(c[l][dummy] + tmp[l]);
	}
      }

      
    }
   
    // printf("A end\n");
#ifdef CE_PARALLEL
    if (err_occurred)
      GERR1("covariance model returns NAs.%.90s", CONTACT);
#endif

    //   int zz=0; printf("\n");for (int i=0; i<2*mtot; i++) for (int j=0; j<4; j++) if (c[j][i]!=0 && ((++zz) < 100000)) printf("%d%d%10g ", i,j,c[j][i]); printf("\n");

    // printf("dim=%d\n", dim);

#ifdef CE_PARALLEL
#define S_FFT(l) s->FFT + l
#define END vdim
#else
#define S_FFT(l) &(s->FFT)
#define END 1    
#endif     
    for (int j=0; j<END; j++) {
      for (int k=0; k<END; k++) {
	//
 	if (k >= j) { // nicht aendern; bei do_ce brauche ich die erste Spalte von s->FFT[] gefuellt 
	  int VARIABLE_IS_NOT_USED l= j*vdim + k; // unused in some cases!!
	  assert(S_FFT(l) != NULL);
	  //	  printf("hier %d %d %d %d\n", j, k, vdim, l);
	  if ((err = fastfourierInit(mm, dim, S_FFT(l))) != NOERROR)
	    goto ErrorHandling;
	  //	  printf("ende\n");
	}
      }
    }

    
    err_occurred = NOERROR;

#ifdef CE_PARALLEL
#pragma omp parallel for num_threads(CORES) if (vdim > 1) collapse(2) schedule(dynamic, 1) reduction(+:err_occurred)
#endif
    for (int j=0; j<vdim; j++) {
      for (int k=0; k<vdim; k++) {
	//
 	if (k >= j) { // nicht aendern; bei do_ce brauche ich die erste Spalte von s->FFT[] gefuellt 
	  int l= j*vdim + k;
	  /*
	  assert({bool ok__ = true;			\
	      for (int ii =0; ii<mtot * 2; ii++) {		\
		if (ISNAN(c[l][ii])) {				\
		  PRINTF("%d %d %10g\n", l, ii, c[l][ii]);	\
		  ok__ = false;					\
		}						\
	      }; ok__;						\
	    });
	  */
#ifdef CE_PARALLEL 
#else	  
	  if (PL>=PL_STRUCTURE) { LPRINT("FFT.0.."); }
#endif
	  err_occurred += fastfourier(c[l], mm, dim, !false, S_FFT(l));

	  //	  for (int t=0; t<mtot; t++) printf("%f ", c[l][t]);printf("\n");

	  
#ifdef CE_PARALLEL
#else	  
   	  if (err_occurred != NOERROR) goto ErrorHandling;
	  if (PL>=PL_STRUCTURE) { LPRINT("done %d.\n", l); }
#endif
	  //	if (false) {
	  //	  for (int ii =0; ii<mtot * 2; ii++) 
	  //	    if (ISNAN(c[l][ii])) GERR("fourier transform returns NAs");
	  //	}
        } 
      }
    }
    
#ifdef CE_PARALLEL
    if (err_occurred != NOERROR) { GERR("error within fft"); }
#endif			       
  //printf("AA end\n");
  
    if (PL>=PL_STRUCTURE) { LPRINT("finished\n"); }

    // here the A(k) have been constructed; A(k) =      c[0][k], c[1][k], ...,  c[p-1][k]
    //                                                  c[p][k], ...            c[2p-1][k]
    //                                                  ...                     ...
    //                                                  c[(p-1)*p][k], ...,     c[p*p-1][k]
    // NOTE: only the upper triangular part of A(k) has actually been filled; A is Hermitian.

    // Now:
    // A(k) is placed into R (cf W&C 1999)
    // application of zheev puts eigenvectors
    // of A(k) into R


    if(Lambda!=NULL) {
      for(int l = 0; l<vdim; l++) FREE(Lambda[l]);
      UNCONDFREE(Lambda);
    }
    if( (Lambda = (double **) MALLOC(vdim * sizeof(double *))) == NULL) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
    }
    for(int l=0; l<vdim; l++) {
      if( (Lambda[l] = (double *) MALLOC(mtot * sizeof(double))) == NULL) {
	err=ERRORMEMORYALLOCATION; goto ErrorHandling;
      }
    }
  
	
    Long notposdef = 0;
    if (vdim==1) { // in case vdim==1, the eigenvalues we need are already contained in c[0]
      // parallelisierung lohnt hier nicht!
      for(int i = 0; i<mtot; i++) {
	Lambda[0][i] = c[0][2*i];
	notposdef += FABS(c[0][2*i+1]) > tolIm || Lambda[0][i] <= tolRe;
	// if (!s->positivedefinite) printf("%10e %10e\n", c[0][2*i+1], tolIm);
	c[0][2*i] = 1.0;
	c[0][2*i+1] = 0.0;
      }
    } else {
      int index1, index2, Sign;      
      //  printf("ABC\n");
#ifdef CE_PARALLEL
#pragma omp parallel for num_threads(CORES) if (vdim > 1L) reduction(+:notposdef)
#endif      
      for(int i = 0; i<mtot; i++) {
      // printf("%d %d ", i, mtot);
	int msg,
	  lwork = worksize,
	  twoi=2*i,
	  twoi_plus1=twoi+1;
	complex optim_lwork;

	// construct R (vdim x vdim matrix) (which is idential to A[i] (cf W&C))
	for (int k=0; k<vdim; k++) {
	  for (int l=0; l<vdim; l++) {   
	    index1 = vdim * k + l;
	    if(l>=k) { 
	      index2=index1;
	      Sign= -1; // 10.2.21
	      notposdef += k==l && FABS(c[index2][twoi_plus1]) > tolIm;
	      //	      c[index2][twoi_plus1].r = 0.0; ?!!!
	    }  // obtain lower triangular bit of c as well with
	    else{
	      index2= vdim * l + k;
	      Sign= +1;// 10.2.21
	    }  // (-1) times imaginary part of upper triangular bit

	    R[index1].r = c[index2][twoi];
	    R[index1].i = Sign*c[index2][twoi_plus1];

	    //  printf("%d:%f+%fi ", i,  R[index1].r, R[index1].i);
	  }
	}

	// diagonalize R via Fortran call zheev
	// initialization
	if(false && i==0) {
#define INI_LWORK -1	  
	  lwork = INI_LWORK; // for initialization call to get optimal size of array 'work'

	  F77_CALL(zheev)("V", "U", &vdim, R, // Eigen,cmplx
			  &vdim, // integer
			  tmpLambda,  // double
			  &optim_lwork, // complex
			  &lwork, rwork, &msg); 



           //// martin 29.12.: nur einmal aufrufen]
	  ////                                ganz weit vorne, oder haengt
	  //// optim_lwork.r  von den Werten der Matrix ab ??

	  lwork = (int) optim_lwork.r;
	  //	  printf("lwork = %d %d\n", lwork, worksize); 
	  if (lwork > worksize) BUG;

	  //  INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
	  //      NB = ILAENV(          1, 'ZHETRD', UPLO, N, -1, -1, -1 )
	  // IC = Z, IZ = z, 
  //     C2 = SUBNAM( 2: 3 )
  //    C3 = SUBNAM( 4: 6 )
	  //      C4 = C3( 2: 3 )
			      
	}

	// optim vdim * 

	F77_CALL(zheev)("V", "U", &vdim, R, &vdim, // Eigen
			tmpLambda, work, &lwork, rwork, &msg);



	for (int l=0;l<vdim;l++) {
	  Lambda[l][i] = tmpLambda[l];
	  // printf("%f ", Lambda[l][i]);
	  notposdef += Lambda[l][i] < tolRe;
	}

	// now R contains the eigenvectors of A[i]
	// and Lambda[l][i], l=0,...,vdim-1, contains the eigenvalues of A[i]

	// the following step writes the vdim vdim-dimensional eigenvectors just computed into c


	for(int j=0, k=0; k<vdim; k++) {
	  int endfor = vdimSq + k;
	  for(int l=k; l<endfor; l+=vdim, j++) {
	    c[j][twoi] = R[l].r;
	    c[j][twoi_plus1] = R[l].i;
	  }
	}

	// peter 28.12. doch wieder
	if(false)
	  for(int l=0; l<vdimSq; l++) {
	    c[l][twoi] = R[l].r;
	    c[l][twoi_plus1] = R[l].i;
	  }


	// next: check if all eigenvalues are >= 0 (non-negative definiteness)
	// and do check this within this loop... in case of some lambda <0 there is no point in going
	// on diagonalizing if there are free trial... or not? Hmm.
      }
    }
    s->positivedefinite = !notposdef;

    // if (s->positivedefinite) {
    //  for (int i=0; i<mtot;i++) for(int l=0;l<vdim;l++)  printf("l=%d i=%d %f\n", l, i, Lambda[l][i]); } else printf("NOT POS DEF !!!\n");

    // printf("EEE\n");
    // check if positive definite. If not: enlarge and restart 
    if (!force || (s->trials<trials)) {

      if (equalsnowVariogram(next)) { // 2.2.19 war isnowVariogram. Warum??
	for (int l=0; l<vdim; l++) {
	  if (Lambda[l][0] < 0.0) Lambda[l][0]=0.0;
	}
      } 

      for (int i=0; i<mtot; i++) { // MULT

	// Check if all eigenvalues of A(i) are > 0:
	int ll=0;
	while( (ll<vdim) && (s->positivedefinite &= Lambda[ll][i] >= tolRe) )  
	  // eigenvalues ARE real: drop former line
	{ll++;}                                                      // s->positivedefinite = (lambda.real>=cepar->tol_re) && (lambda.imag<=cepar->tol_im)


	if ( !s->positivedefinite) {
	  int *sum = NULL;
	  double *smallest=NULL;
	  Long *idx = NULL;
	  if (PL>=PL_BRANCHING) {
	    if (Lambda[ll][i] < 0.0) {
	      LPRINT("There are complex eigenvalues\n");
	    } else if (vdim==1) {
	      LPRINT("non-positive eigenvalue: Lambda[%d]=%10e (< %10e).\n",
		     i, Lambda[ll][i], tolRe);
	    } else {
	      LPRINT("non-pos. eigenvalue in dim %d: Lambda[%d][%d]=%10e.\n",
		     ll, i, ll, Lambda[ll][i]);
	    }
	  }
	  if (PL>=PL_ERRORS) { // just for printing the smallest 
	      //                            eigenvalue (min(c))
 	    if( (sum = (int *) CALLOC(vdim, sizeof(int))) == NULL ||
		(smallest = (double *)CALLOC(vdim, sizeof(double))) == NULL || 
		(idx = (Long *) CALLOC(vdim, sizeof(Long))) == NULL) {
	      FREE(sum);
	      FREE(smallest);
	      FREE(idx);
	      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
	    }
	    char percent[]="%";
	    for (int j=i; j<mtot; j++) {
	      for (int k=0; k<vdim; k++) {
		if(Lambda[k][j] < 0) {
		  sum[k]++;
		  if(Lambda[k][j] < smallest[k]) {
		    idx[k] = j;
		    smallest[k] = Lambda[k][j];}
		}
	      }
	    }
	    if(vdim==1) {
	      LPRINT("   There are %d negative eigenvalues (%5.3f%s).\n",
		  sum[0], 100.0 * (double) sum[0] / (double) mtot, percent);
	      LPRINT("   The smallest eigenvalue is Lambda[%ld] =  %10e\n", idx[0], smallest[0]);
	    }
	    else {
	      for (int l=0; l<vdim; l++) {
		LPRINT("   There are %d negative eigenvalues in dimension %d (%5.3f%s).\n",
		    sum[l], l, 100.0 * (double) sum[l] / (double) mtot, percent);
		LPRINT("   The smallest eigenvalue in dimension %d is Lambda[%ld][%d] =  %10e\n", l, idx[l], l, smallest[l]);
	      }
	    }
	    FREE(sum);
	    FREE(smallest);
	    FREE(idx);
	  }
	  break; // break loop 'for(i=0 ...)'
	}	
      }

      //print("%d %d %d %d\n", s->positivedefinite, s->trials, trials, force);
     
      if (!s->positivedefinite && (s->trials<trials)) { 
#ifdef CE_PARALLEL
	for (int i=0; i<vdim * vdim; i++) FFT_destruct(s->FFT + i);
#else
	FFT_destruct(&(s->FFT));
#endif	        
	switch (strategy) {
	  case 0 :
	    for (int i=0; i<dim; i++) { /* enlarge uniformly in each direction, maybe 
				       this should be modified if the grid has 
				       different sizes in different directions */
		hilfsm[i] *= multivariate ? 3 : 2;
	    }
	    break;
	  case 1 :  
	    double cc, maxcc, hx[MAXCEDIM];

	    int maxi;
	    maxi = UNSET;
	    maxcc = RF_NEGINF;
	    for (int i=0; i<dim; i++) hx[i] = 0.0;
	    for (int i=0; i<dim; i++) {

	      // MUSS HIER ETWAS GEAENDERT WERDEN FUER MULTIVARIAT?

	      hx[i] = steps[i] * mm[i]; 

	      cc=0;
	      COV(hx, info, next, tmp);
	      for (int l=0; l<vdimSq; l++) cc += FABS(tmp[l]);
	      //if (PL>2) { LPRINT("%d cc=%10e (%10e)",i,cc,hx[i]); }
	      if (cc>maxcc) {
		maxcc = cc;
		maxi = i; 
	      }
	      hx[i] = 0.0;
	    }
	    assert(maxi>=0); 
	    hilfsm[maxi] *= multivariate ? 3 : 2;
	    break;
	  default:
	    GERR("unknown strategy for circulant embedding");
	}
      }
    } else {if (PL>=PL_BRANCHING) {LPRINT("forced.\n");} }
    R_CheckUserInterrupt();

    //   printf("trials=%d %d %d\n", s->trials, trials, s->positivedefinite);

  } // while (!s->positivedefinite && (s->trials<trials)) 354-706


  assert(mtot>0);
  if (s->positivedefinite || force) { 
    // correct theoretically impossible values, that are still within 
    // tolerance CIRCEMBED.tol_re/CIRCEMBED.tol_im 

    double r;
    r = Lambda[0][0]; // MULT

    // keine Parallelisierung!

    for (int i=0; i<mtot;i++) {
      int
	twoi= 2 * i,
	twoi_plus1 = twoi+1;
      for(int l=0;l<vdim;l++) {
	if(Lambda[l][i] > 0.0) Lambda[l][i] = SQRT(Lambda[l][i]); 
	else {
	  if(Lambda[l][i] < r) r = Lambda[l][i];
	  Lambda[l][i] = 0;
	}
      }

      // now compute R[i]*Lambda^(1/2)[i]
      for(int j=0;j<vdim;j++) {
	for(int l=0; l<vdim; l++) {
	  c[j*vdim+l][twoi] *= Lambda[l][i];
	  c[j*vdim+l][twoi_plus1] *= Lambda[l][i];
	}
      }
    }
    if (PL>=PL_BRANCHING) {
      if (r < -GENERAL_PRECISION) {
	LPRINT("using approximating circulant embedding:\n");
	LPRINT("\tsmallest real part has been %10e \n", r);
      }
    }

  } else {
    //printf("FEHLER\n");
    GERR("embedded matrix has (repeatedly) negative eigenvalues and approximation is not allowed (force=FALSE)");
  }
  if (PL >= PL_SUBDETAILS) {
    for (int i=0; i < 2 * mtot; i++)
      for (int l=0; l<vdimSq; l++)
	{LPRINT("%10g ",c[l][i]);} 
    LPRINT("\n");
  }


  s->dependent = dependent;
  s->new_simulation_next = true;
  for(int i=0; i<dim; i++) { // set anyway -- does not cost too much
    s->cur_square[i] = 0;
    s->max_squares[i] = hilfsm[i] / s->nn[i];
    s->square_seg[i] = cumm[i] * 
      (s->nn[i] + (hilfsm[i] - s->max_squares[i] * s->nn[i]) / 
       s->max_squares[i]);
  }

  if ( (s->d=(double **) CALLOC(vdim, sizeof(double *))) == NULL ) {
      err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  for(int l=0; l<vdim;l++) {
      if( (s->d[l] = (double *) CALLOC(2L * mtot, sizeof(double))) == NULL) {
	  err=ERRORMEMORYALLOCATION;goto ErrorHandling;
      }
  }

  if ((s->gauss1 = (complex *) MALLOC(sizeof(complex) * vdim)) == NULL) {
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }
  if ((s->gauss2 = (complex *) MALLOC(sizeof(complex) * vdim)) == NULL) {
    err=ERRORMEMORYALLOCATION; goto ErrorHandling;
  }

  err = ReturnOwnField(cov);

ErrorHandling:
  s->c = c;
  if (Lambda != NULL) {
    for(int l = 0; l<vdim; l++) FREE(Lambda[l]);
    UNCONDFREE(Lambda);
  }

  FREE(tmp);
  FREE(tmpLambda);
  FREE(rwork);
  FREE(R);
  FREE(work);
  

//  printf("end init\n"); BUG;
  cov->simu.active = err == NOERROR;
  RETURN_ERR(err);
}


void kappa_ce(int i, model *cov, int *nr, int *nc){
  kappaGProc(i, cov, nr, nc);
  if (i == CE_MMIN) *nr = 0;
}


int check_ce_basic(model *cov) { 
  globalparam *global = &(cov->base->global);
  // auf co/stein ebene !
  //  model *next=cov->sub[0];
  //  location_type *loc = Loc(cov);
  int i, //err,
    dim = ANYDIM; // taken[MAX DIM],
  ce_param *gp  = &(global->ce); // ok
  
  FRAME_ASSERT_GAUSS_INTERFACE;
  ASSERT_CARTESIAN;
  ASSERT_UNREDUCED;

  kdefault(cov, CE_FORCE, (int) gp->force);
  if (PisNULL(CE_MMIN)) {
    PALLOC(CE_MMIN, dim, 1);
    for (i=0; i<dim; i++) {
      P(CE_MMIN)[i] = gp->mmin[i];
    }
  } 
  kdefault(cov, CE_STRATEGY, (int) gp->strategy);
  kdefault(cov, CE_MAXGB, gp->maxGB);
  kdefault(cov, CE_MAXMEM, (int) gp->maxmem);
  kdefault(cov, CE_TOLIM, gp->tol_im);
  kdefault(cov, CE_TOLRE, gp->tol_re);
  kdefault(cov, CE_TRIALS, gp->trials);
  kdefault(cov, CE_USEPRIMES, (int) gp->useprimes);
  kdefault(cov, CE_DEPENDENT, (int) gp->dependent);
  kdefault(cov, CE_APPROXSTEP, gp->approx_grid_step);
  kdefault(cov, CE_APPROXMAXGRID, gp->maxgridsize);

  RETURN_NOERROR;
}

int check_ce(model *cov) { 
  model *next=cov->sub[0];
  //  location_type *loc = Loc(cov);
  int err,
    dim = ANYDIM;

  //  PMI0(cov); printf("%d %s\n", hasInterfaceFrame(cov->calling), NAME(cov->calling));
  FRAME_ASSERT_GAUSS_INTERFACE; 
  ASSERT_UNREDUCED;
  ASSERT_ONESYSTEM;
  RESERVE_BOXCOX;
  
  if (dim > MAXCEDIM) RETURN_ERR(ERRORCEDIM);
  if ((err = check_ce_basic(cov)) != NOERROR) RETURN_ERR(err);
  if ((err = checkkappas(cov, false)) != NOERROR) RETURN_ERR(err);
  if (Loctsdim(cov) > MAXCEDIM || OWNTOTALXDIM > MAXCEDIM)
    RETURN_ERR(ERRORCEDIM);
  
  if (cov->key != NULL) {
    if ((err = CHECK_PASSFRAME(cov->key, GaussMethodType)) != NOERROR) {
       //PMI(cov->key); XERR(err);
       //printf("specific ok\n");
       // crash();
       RETURN_ERR(err);
     }
  } else {
    if ((err = CHECK(next, dim, dim, PosDefType, XONLY, CARTESIAN_COORD,
		     SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
      //      APMI(cov);
      //    XERR(err); APMI(cov);
      if ((err = CHECK(next, dim, dim, VariogramType, XONLY, SYMMETRIC, 
		       SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
	RETURN_ERR(err);
      }
    }
    if (next->pref[CircEmbed] == PREF_NONE) {
      //PMI(cov);
      //    assert(false);
      RETURN_ERR(ERRORPREFNONE);
    }
    if (!isPosDef(NEXTTYPE(0))) SERR("only covariance functions allowed.");
  }
  setbackward(cov, next);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
   if ((err = checkkappas(cov, true)) != NOERROR) RETURN_ERR(err);
  RETURN_NOERROR;
}


void range_ce(model VARIABLE_IS_NOT_USED *cov, int k, int i, int j,
	      simple_range_type *range) {
  switch(k) {
    GAUSS_BOXCOX_RANGE;
    BooleanRange(CE_FORCE);
  case CE_MMIN : 
    range->min = RF_NEGINF; 
    range->max = RF_INF;
    range->pmin = -3;
    range->pmax = 10000;
    range->openmin = true;
    range->openmax = true;
    break;

    BooleanRange(CE_STRATEGY);

  case CE_MAXGB :
    range->min = 0; 
    range->max = RF_INF;
    range->pmin = 0;
    range->pmax = 10;
    range->openmin = false;
    range->openmax = false;
    break;
  case CE_MAXMEM :
    range->min = 0; 
    range->max = MAXINT;
    range->pmin = 0;
    range->pmax = 50000000;
    range->openmin = false;
    range->openmax = false;
    break;
  case CE_TOLIM :
    range->min = 0.0; 
    range->max = RF_INF;
    range->pmin = 0.0;
    range->pmax = 1e-2;
    range->openmin = false;
    range->openmax = false;
    break;
  case CE_TOLRE :
    range->min = RF_NEGINF; 
    range->max = RF_INF;
    range->pmin = 1e-2;
    range->pmax = 0;
    range->openmin = false;
    range->openmax = false;
    break;
  case CE_TRIALS :
    range->min = 1; 
    range->max = 10;
    range->pmin = 1;
    range->pmax = 3;
    range->openmin = false;
    range->openmax = true;
    break;
    
    BooleanRange(CE_USEPRIMES);
    BooleanRange(CE_DEPENDENT);
  
  case CE_APPROXSTEP :
    range->min = 0.0; 
    range->max = RF_INF;
    range->pmin = 1e-4;
    range->pmax = 1e4;
    range->openmin = true;
    range->openmax = true;
    
    break;
  case CE_APPROXMAXGRID :
    range->min = 1; 
    range->max = RF_INF;
    range->pmin = 100;
    range->pmax = 5e7;
    range->openmin = false;
    range->openmax = false;
    break;
  default : BUG;
  }
}

 

void do_circ_embed(model *cov, gen_storage VARIABLE_IS_NOT_USED *S){
  globalparam *global = &(cov->base->global);
  //printf("start do\n");
  assert(cov != NULL);
  assert(Locgrid(cov));
  SAVE_GAUSS_TRAFO;
      
  int  HalfMp1[MAXCEDIM], HalfMaM[2][MAXCEDIM], index[MAXCEDIM], 
    dim, *mm, *cumm, *halfm, // MULT added vars l, pos
    vdim;
  double invsqrtmtot,
    **c=NULL,
    **d=NULL,
    *dd=NULL,
    *res = cov->rf; // MULT *c->**c
  bool vfree[MAXCEDIM+1], noexception; // MULT varname free->vfree
  Long mtot, start[MAXCEDIM], end[MAXCEDIM];
getStorage(s ,   ce); 
  location_type *loc = Loc(cov);
  
  complex *gauss1 = s->gauss1, *gauss2 = s->gauss2;

  if (s->stop) XERR(ERRORNOTINITIALIZED);
 
  /* 
     implemented here only for rotationsinvariant covariance functions
     for arbitrary dimensions;
     (so it works only for even covariance functions in the sense of 
     Wood and Chan,p. 415, although they have suggested a more general 
     algorithm;) 
     
  */  
  
  dim = OWNLOGDIM(0);
  mm = s->m;
  cumm = s->cumm;
  halfm = s->halfm;
  c = s->c;
  d = s->d;
  mtot= s->mtot;
  vdim = s->vdim; // so p=vdim
  
#ifdef CE_PARALLEL
  omp_set_num_threads(CORES);
#endif
  
  for (int i=0; i<dim; i++) {
    HalfMp1[i] =  ((mm[i] % 2)==1) ? -1 : halfm[i];
    HalfMaM[0][i] = halfm[i];
    HalfMaM[1][i] = mm[i] - 1;
  }


  //if there are doubts that the algorithm below reaches all elements 
  //of the matrix, uncomment the lines containing the variable xx
  //
  //bool *xx; xx = (bool*) MALLOC(sizeof(bool) * mtot);
  //for (i=0; i<mtot;i++) xx[i]=true;
  
  invsqrtmtot = 1.0 / SQRT((double) mtot);

  if (PL>=PL_STRUCTURE) { LPRINT("Creating Gaussian variables... \n"); }
  /* now the Gaussian r.v. have to defined & multiplied with SQRT(FFT(c))*/

  for (int i=0; i<dim; i++) {
    index[i]=0; 
    vfree[i] = false;
  }
  vfree[MAXCEDIM] = false;

  // PM 12/12/2008
  // don't simulate newly in case this is an even call
  s->cur_call_odd = !(s->cur_call_odd);
  s->new_simulation_next = (s->new_simulation_next && s->cur_call_odd);

  if (s->new_simulation_next) {

    // MULT
 
    
    // martin: 13.12. check
    //for(k=0; k<vdim; k++) {
    //      dd = d[k];
    //  for (j=0; j<mtot; j++) dd[j] = RF_NA;
    //    }

    int i, j;
    for(;;) {      
      i = j = 0;
      noexception = false;
      for (int k=0; k<dim; k++) {
	i += cumm[k] * index[k];
	if ((index[k]==0) || (index[k]==HalfMp1[k])) {
	  j += cumm[k] * index[k];
	} else {
	  noexception = true; // not case 2 in prop 3 of W&C
	  j += cumm[k] * (mm[k] - index[k]);
	}
      }
      
      //A 0 0 0 0 0 0.000000 0.000000
      //B 0 0 0 0 0 -26.340334 0.000000

      //    if (index[1] ==0) print("A %d %d %d %d %d %10g %10g\n", 
      //	   index[0], index[1], i, j, noexception, d[i], d[i+1]);
      //      if (PL>=10) { LPRINT("cumm..."); }
      i <<= 1; // since we have to index imaginary numbers
      j <<= 1;

      //      print("ij %d %d (%d)\n", i, j, noexception);

      if (noexception) { // case 3 in prop 3 of W&C

        for(int k=0; k<vdim; k++) {
	  gauss1[k].r = GAUSS_RANDOM(1.0);
	  gauss1[k].i = GAUSS_RANDOM(1.0);
	  gauss2[k].r = GAUSS_RANDOM(1.0);
	  gauss2[k].i = GAUSS_RANDOM(1.0);
       }

	
        for(int k=0; k<vdim; k++) {
	  dd = d[k];
	  dd[j] = dd[i] = dd[j+1] = dd[i+1] = 0.0;// martin:13.12. neu eingefuegt
	  for(int l=0; l<vdim; l++) {
	    int pos = k * vdim + l;
	    dd[i] += c[pos][i] * gauss1[l].r - c[pos][i+1] * gauss1[l].i;   
	    // real part
	    dd[i+1] += c[pos][i+1] * gauss1[l].r + c[pos][i] * gauss1[l].i; 
	    // imaginary part 
	    
	    dd[j] += c[pos][j] * gauss2[l].r - c[pos][j+1] * gauss2[l].i;  
	    // real part    // martin: 13.12. geaendert von = + zu +=
	    dd[j+1] += c[pos][j+1] * gauss2[l].r + c[pos][j] * gauss2[l].i;  
	    // imaginary part
	  }
        }

      } else { // case 2 in prop 3 of W&C

        for(int k=0; k<vdim; k++) {
	  gauss1[k].r = GAUSS_RANDOM(1.0);
	  gauss1[k].i = GAUSS_RANDOM(1.0);
	}

        for(int k=0; k<vdim; k++) {
	  dd = d[k];
	  dd[i+1] = dd[i] = 0.0; //  martin: 13.12. neu eingefuegt
	  for(int l=0; l<vdim; l++) {
	    int pos = k * vdim + l;
	    dd[i] += c[pos][i] * gauss1[l].r - c[pos][i+1] * gauss1[l].i;   
	    // real part
	    dd[i+1] += c[pos][i+1] * gauss1[l].r + c[pos][i] * gauss1[l].i;  
	    // imaginary part
	  }
        }	
      }

      /* 
	 this is the difficult part.
	 We have to run over roughly half the points, but we should
	 not run over variables twice (time lost)
	 Due to case 2, we must include halfm.
	 
	 idea is:
	 for (i1=0 to halfm[dim-1])
	 if (i1==0) or (i1==halfm[dim-1]) then endfor2=halfm[dim-2] 
	 else endfor2=mm[dim-2]
         for (i2=0 to endfor2)
	 if ((i1==0) or (i1==halfm[dim-1])) and
	 ((i2==0) or (i2==halfm[dim-2]))
	 then endfor3=halfm[dim-3] 
	 else endfor3=mm[dim-3]
	 for (i3=0 to endfor3)
	 ....
	 
	 i.e. the first one that is not 0 or halfm (regarded from dim-1 to 0)
	 runs over 0..halfm, all the others over 0..m
	 
	 this is realised in the following allowing for arbitrary value of dim
	 
	 vfree==true   <=>   endfor==mm[]	 
      */
      int k=0; 
      if (++index[k]>HalfMaM[vfree[k]][k]) {
	// in case k increases the number of indices that run over 0..m increases
	vfree[k] = true;
	index[k]= 0; 
	k++;
	while((k<dim) && (++index[k]>HalfMaM[vfree[k]][k])) {
	  vfree[k] = true;
	  index[k]= 0; 
	  k++;
	}
	if (k>=dim) break;
	// except the very last (new) number is halfm and the next index is
	// restricted to 0..halfm
	// then k decreases as Long as the index[k] is 0 or halfm
	if (!vfree[k] && (index[k]==halfm[k])){//index restricted to 0..halfm?
	  // first: index[k] is halfm? (test on ==0 is superfluent)
	  k--;
	  while ( (k>=0) && ((index[k]==0) || (index[k]==halfm[k]))) {
	    // second and following: index[k] is 0 or halfm?
	    vfree[k] = false;
	    k--;
	  }
	}
      }
    }
      
    // MULT
    //
#ifdef CE_PARALLEL
#pragma omp parallel for num_threads(CORES) schedule(static,1)
#endif	    
    for(int k=0; k<vdim; k++) {
      //      printf("ABCD %d %d\n", vdim, omp_get_num_threads()); 
     //
      //printf("k=%d %10g\n", k);
      fastfourier(d[k], mm, dim, !true, S_FFT(k));
      //printf("end do");
    }
  } // if ** simulate **
  // printf("ABCD end\n"); 



  /* now we correct the result of the fastfourier transformation
     by the factor 1/SQRT(mtot) and read the relevant matrix out of 
     the large vector c */
  for(int i=0; i<dim; i++) {
    index[i] = start[i] = s->cur_square[i] * s->square_seg[i];
    end[i] = start[i] + cumm[i] * s->nn[i];
    //    print("%d start=%d end=%d cur=%d max=%d seq=%d cumm=%d, nn=%d\n",
    //	   i, (int) start[i], (int) end[i],
    //	   s->cur_square[i], s->max_squares[i],
    //	   (int) s->square_seg[i], cumm[i], s->nn[i]);
  }
  Long totpts = loc->totalpoints; // MULT

  /*
    for(l=0; l<vdim; l++) { // MULT
    for(k=0; k<dim; k++) index[k] = start[k];

    for (i=l*totpts; i<(l+1)*totpts; i++){
    j=0; for (k=0; k<dim; k++) {j+=index[k];}
    //res[i] += d[l][2 * j + !(s->cur_call_odd)] * invsqrtmtot;
    res[i] += d[l][2 * j] * invsqrtmtot;
    for(k=0; (k<dim) && ((index[k] += cumm[k]) >= end[k]); k++) {
    index[k]=start[k];
    }
    }
    }
  */
  // 


  for(int k=0; k<dim; k++) index[k] = start[k];

  bool vdim_close_together = global->general.vdim_close_together;
  for (int L=0; L<totpts; L++) {
    int j=0;
    for (int k=0; k<dim; k++) {j+=index[k];}
    int xvdim, xtotpts;
    if (vdim_close_together) {
      xvdim = vdim;
      xtotpts = 1;
    } else {
      xvdim = 1;
      xtotpts = totpts;
    }
      
    for (int l=0; l<vdim; l++) { // MULT
      res[xvdim * L + l * xtotpts] =(double)(d[l][2 * j + !(s->cur_call_odd)] * invsqrtmtot);
    }
    for(int k=0; (k<dim) && ((index[k] += cumm[k]) >= end[k]); k++) {
      index[k]=start[k];      
    }
  }

 

  s->new_simulation_next = true;
  if (s->dependent && !(s->cur_call_odd)) {
    int k=0; 
    while(k<dim && (++(s->cur_square[k]) >= s->max_squares[k])) {
      s->cur_square[k++]=0;
    }    
    s->new_simulation_next = k == dim; 
  }
  //
  //  print("dep=%d odd=%d squ=%d %d\n", s->dependent, s->cur_call_odd,
  //	 s->cur_square[0],s->cur_square[1]);
  s->stop |= s->new_simulation_next && s->d == NULL;

  BOXCOX_INVERSE;
  
}



int GetOrthogonalUnitExtensions(double * aniso, int dim, double *grid_ext) {
  int k,i,j,l,m, job=01, Err, dimsq, ev0, jump, endfor;
  double *s=NULL, G[MAXCEDIM+1], e[MAXCEDIM], D[MAXCEDIM], *V=NULL;
  dimsq = dim * dim;

  assert(aniso != NULL);
  
  s = (double*) MALLOC(dimsq * sizeof(double));
  V = (double*) MALLOC(dimsq * sizeof(double));
  for (k=0; k<dim; k++) {
    // kte-Zeile von aniso wird auf null gesetzt/gestrichen -> aniso_k
    // s = aniso %*% aniso_k 
    // s hat somit mindestens Eigenwert der 0 ist.
    //der zugeheorige EVektor wird miit der k-ten Spalte von aniso multipliziert
    // und ergibt dann den Korrekturfaktor
    for (i=0; i<dim; i++) {
      for (l=j=0; j<dimsq; j+=dim) {
	s[j + i] = 0.0;
	jump = l + k;
	endfor = l + dim;
	for (m=i; l<endfor; l++, m += dim) {
	  if (l!=jump) s[i + j] +=  aniso[m] * aniso[l];
	}
      }
    }
    //    for (print("\n\n s\n"), j=0; j<dim; j++, print("\n")) 
    //	for (i=0; i<dimsq; i+=dim) print("%10g ", s[i+j]);
    //    for (print("\n\n aniso\n"), j=0; j<dim; j++, print("\n")) 
    //	for (i=0; i<dimsq; i+=dim) print("%10g ", aniso[i+j]);	
    
    // note! s will be distroyed by dsvdc!
    F77_NAME(dsvdc)(s, &dim, &dim, &dim, D, // SVD
		     e, NULL /* U */, &dim, V, &dim, G, &job, &Err);
    if (Err != NOERROR) { Err = XERRORSVD; goto ErrorHandling;}
    ev0 = UNSET;
    for (i=0; i<dim; i++) {
      if (FABS(D[i]) <= EIGENVALUE_EPS) {
	if (ev0 == UNSET) ev0=i;
	else {Err=XERRORANISO_INV; goto ErrorHandling;}
      }
    }
    
    //    for (print("\n\n V\n"), j=0; j<dim; j++, print("\n")) 
    //	for (i=0; i<dimsq; i+=dim) print("%10g ", V[i+j]);

    grid_ext[k] = 0.0;
    ev0 *= dim;
    for (i=0; i<dim; i++) {
      //	print("%d : %10g %10g %10g\n", k,
      //	       V[ev0 + i], aniso[k + i * dim], V[ev0 + i] / aniso[k + i * dim]);
      grid_ext[k] += V[ev0 + i] * aniso[k + i * dim];
    }
    grid_ext[k] = FABS(grid_ext[k]);
  }
  FREE(V);
  FREE(s);
  return NOERROR;

 ErrorHandling:
  FREE(V);
  FREE(s);
  return Err;
  
}






int check_local_proc(model *cov) {
  int 
    dim = ANYDIM,
    err=NOERROR;
  model 
    *key = cov->key,
    *next = cov->sub[0],
    *sub = key != NULL ? key : next;
  //  location_type *loc = Loc(cov);
  bool cutoff = COVNR == CE_CUTOFFPROC_USER || COVNR==CE_CUTOFFPROC_INTERN;
  if (!cutoff && COVNR!=CE_INTRINPROC_USER && COVNR!=CE_INTRINPROC_INTERN)
    BUG;

  FRAME_ASSERT_GAUSS_INTERFACE;
  ASSERT_UNREDUCED;
  ASSERT_ONESYSTEM;
  RESERVE_BOXCOX;
  
   if ((err = check_ce_basic(cov)) != NOERROR) RETURN_ERR(err);
   if (dim > MAXCEDIM) RETURN_ERR(ERRORCEDIM);
 
  if (key != NULL) {
    // falls nicht intern muessen die parameter runter und rauf kopiert
    // werden, ansonsten sind die kdefaults leer.
    // falls hier gelandet, wird RP*INTERN von RP* aufgerufen.
    // oder RP*INTERN wird Circulant aufrufen
    model *intern = cov,
      *RMintrinsic = key->sub[0];
    while (intern != NULL && MODELNR(intern) != CE_INTRINPROC_INTERN &&
	   MODELNR(intern) != CE_CUTOFFPROC_INTERN) {
      // normalerweise ->key, aber bei $ ist es ->sub[0]
      intern = intern->key == NULL ? intern->sub[0] : intern->key;
    }
    if (intern == NULL) {
      BUG;
    }
    else if (intern != cov) // cov not intern, ie user
      paramcpy(intern, cov, true, true, false, false, false);  // i.e. INTERN
    else if (MODELNR(key) == CE_INTRINPROC_INTERN || 
	     MODELNR(key) == CE_CUTOFFPROC_INTERN) { 
      // 2x intern hintereinander, d.h. es wird eine approximation durchgefuehrt
      paramcpy(key, cov, true, true, false, false, false);        
    } else {
      //if (RMintrinsic->nr == GAUSSPROC) RMintrinsic = RMintrinsic->sub[0];
      if (MODELNR(RMintrinsic) != CUTOFF && MODELNR(RMintrinsic) != STEIN) {
	BUG;
      }

      if (!PisNULL(LOCPROC_DIAM)) 
	kdefault(RMintrinsic, pLOC_DIAM, P0(LOCPROC_DIAM));
      if (!PisNULL(LOCPROC_R)) 
	kdefault(RMintrinsic, pLOC_DIAM, P0(LOCPROC_R));
    }
     
    if ((err = CHECK(sub, dim,  dim, ProcessType, KERNEL, CARTESIAN_COORD, 
		     SUBMODEL_DEP, GaussMethodType)) != NOERROR) {
      RETURN_ERR(err);
    }
    if (intern == cov) { // all the other parameters are not needed on the 
      //                    upper levels
      if (PisNULL(LOCPROC_DIAM))
	kdefault(cov, LOCPROC_DIAM, PARAM0(RMintrinsic, pLOC_DIAM));  
    }
  } else {
    if ((err = CHECK(sub, dim,  1, VariogramType,
             //cutoff ? PosDefType : VariogramType,
              XONLY, ISOTROPIC, SUBMODEL_DEP, GaussMethodType))
     != NOERROR) {
       if (isDollar(next) && !PARAMisNULL(next, DANISO)) {
     // if aniso is given then xdimprev 1 does not make sense
     err = CHECK(sub, dim, dim, VariogramType,
             //cutoff ? PosDefType : VariogramType,
             XONLY, ISOTROPIC, SUBMODEL_DEP, GaussMethodType);
       }
       //
       if (err != NOERROR) RETURN_ERR(err);
       // PMI(cov, "out");
       //assert(false);
     }
   }

  //printf("check intern end\n");

  // no setbackward ?!
  setbackward(cov, sub); 
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) != NOERROR) RETURN_ERR(err);
  VDIM0 = VDIM1 = sub->vdim[0];
  
  RETURN_NOERROR;
}

int init_circ_embed_local(model *cov, gen_storage *S){
 // being here, leadin $ have already been put away
  // so a new $ is included below

  model //*dummy = NULL,
    *key = cov->key;
  location_type *loc = Loc(cov);
  int instance, i, d, 
    timespacedim = Loctsdim(cov),
    cncol = ANYDIM,
    err = NOERROR;
  //  localCE_storage *s=NULL;
  double grid_ext[MAXCEDIM], old_mmin[MAXCEDIM];
 
  int first_instance;
  double
    *mmin = P(CE_MMIN);

  double sqrt2a2;
  bool is_cutoff = COVNR == CE_CUTOFFPROC_INTERN;

  ASSERT_ONESYSTEM;

  cov->method = is_cutoff ? CircEmbedCutoff : CircEmbedIntrinsic;

  if (loc->distances) RETURN_ERR(ERRORFAILED);


  if (loc->caniso != NULL) {
    if (loc->cani_ncol != loc->cani_nrow || 
	loc->cani_ncol != timespacedim) RETURN_ERR(ERRORCEDIM);
    err = GetOrthogonalUnitExtensions(loc->caniso, timespacedim, grid_ext);
    if (err != NOERROR) RETURN_ERR(err);
  } else {
    for (i=0; i<timespacedim; i++) grid_ext[i] = 1.0;
  }

  //printf("diameter %d \n",LOCPROC_DIAM); APMI(cov);

  if (!PARAMisNULL(key->sub[0], LOCPROC_R))
    kdefault(key, pLOC_DIAM, P0(LOCPROC_R));
 
  kdefault(key, CE_FORCE, P0INT(CE_FORCE));

  PARAMFREE(key, CE_MMIN);
  PARAMALLOC(key, CE_MMIN, ANYDIM, 1);
  PCOPY(key, cov, CE_MMIN);
  
  kdefault(key, CE_STRATEGY, P0INT(CE_STRATEGY));
  kdefault(key, CE_MAXGB, P0(CE_MAXGB));
  kdefault(key, CE_MAXMEM, P0INT(CE_MAXMEM));
  kdefault(key, CE_TOLIM, P0(CE_TOLIM));
  kdefault(key, CE_TOLRE, P0(CE_TOLRE));
  kdefault(key, CE_TRIALS, P0INT(CE_TRIALS));
  kdefault(key, CE_USEPRIMES, P0INT(CE_USEPRIMES));
  kdefault(key, CE_DEPENDENT, P0INT(CE_DEPENDENT));
  RESERVE_BOXCOX;

  //  APMI(key);
 
  assert(VDIM0 == VDIM1);
  err = CHECK_PASSTF(key, GaussMethodType, VDIM0, GaussMethodType);
  //  err = CHECK(key, cov->tsdim, cov->xdimprev,
  //	      GaussMethodType,
  //            OWNDOM(0), OWNISO(0), SUBMODEL_DEP, GaussMethodType); 
  if ((err < MSGLOCAL_OK && err != NOERROR) 
      || err >=MSGLOCAL_ENDOFLIST // 30.5.13 : neu
      ) RETURN_ERR(err);

  first_instance = err != NOERROR;
  for (d=0; d<timespacedim; d++) old_mmin[d] = mmin[d];
  

  model *local ;
  local = key->sub[0];
  GETSTORAGE(ss , local,   localCE); 
  assert(ss != NULL);
  localvariab *q;
  q = ss->q;
  assert(q != NULL);
  assert(err == NOERROR);
  for (instance = first_instance; instance < 2; instance++) {
     if (instance == 1) {
      if (q->msg != MSGLOCAL_OK) {
	if (!DefList[MODELNR(local)].alternative(local)) break;
      } else {
	if (first_instance == 0 && err != NOERROR) goto ErrorHandling;
	else BUG;
      }
    } else assert(instance == 0);
    
    for (d=0; d<timespacedim; d++) {
      if (old_mmin[d]==0.0) {
	
	mmin[d] = - q->R / 
	  (grid_ext[d] * (loc->xgr[d][XLENGTH] - 1.0) * loc->xgr[d][XSTEP]); // OK
	if (mmin[d] > -1.0) mmin[d] = - 1.0;
	
      }
    }
    if ((err = init_circ_embed(key, S)) == NOERROR) break;
  }
  
  if (err != NOERROR) goto ErrorHandling;

  // printf("\n\n\n\n\n\n\nCOVNR =%d %d %d\n", COVNR, CE_INTRINPROC_INTERN, cov->zaehler);
  if (COVNR == CE_INTRINPROC_INTERN) {
    sqrt2a2 = SQRT(2.0 * q->intrinsic.A2);//s. Stein (2002) timespacedim * cncol
    if (loc->caniso == NULL) ss->correction[0] = sqrt2a2;
    else {
      int rowcol = cncol * cncol;   
      for (i=0; i<rowcol; i++) ss->correction[i] = sqrt2a2 * loc->caniso[i];
    }
  } else {
    assert(COVNR == CE_CUTOFFPROC_INTERN);
  }

  assert(err == NOERROR);
  
  ReturnOtherField(cov, cov->key);
  if ((err = kappaBoxCoxParam(cov, GAUSS_BOXCOX)) !=NOERROR) goto ErrorHandling;
 
 ErrorHandling :
  for (d=0; d<timespacedim; d++) mmin[d] = old_mmin[d];
  cov->simu.active = err == NOERROR;

  //  printf("%d\n", local->SlocalCE != NULL);
  // APMI(cov);
 
  RETURN_ERR(err);
}


int struct_ce_local(model *cov, model VARIABLE_IS_NOT_USED **newmodel) {
  model 
    *next = cov->sub[0];
  int err;
  bool cutoff = COVNR ==  CE_CUTOFFPROC_INTERN;
  
  if (next->pref[cutoff ? CircEmbedCutoff : CircEmbedIntrinsic] == PREF_NONE)
    RETURN_ERR(ERRORPREFNONE);

  if (cov->key != NULL) COV_DELETE(&(cov->key), cov);
  if ((err = covcpy(&(cov->key), next)) != NOERROR) RETURN_ERR(err);
  addModelKey(cov, cutoff ? CUTOFF : STEIN);
  addModelKey(cov, CIRCEMBED);
  
  // crash(cov);

  // da durch init die relevanten dinge gesetzt werden, sollte
  // erst danach bzw. innerhalb von init der check zum ersten Mal
  //APMI(cov);
  RETURN_NOERROR;
}


void do_circ_embed_cutoff(model *cov, gen_storage *S) {  
  globalparam *global = &(cov->base->global);
  double   *res = cov->rf;
  model *key = cov->key,
     *sub = key;
   assert(key != NULL);
   sub = sub->key != NULL ? sub->key : sub->sub[0];
   assert(sub != NULL);
   
   int // k, row = loc->timespacedim,
     vdim = VDIM0;
   Long totpts = Loctotalpoints(cov);


    // printf("--- \n \n do_circ_embed_cutoff \n \n  ---");
    //PMI(sub);
   GETSTORAGE(s , sub,     localCE); 
    assert(s != NULL);
    localvariab *q = s->q;
    assert(q != NULL);


    // PMI(sub);

    do_circ_embed(key, S);

    if ( VDIM0 > 1) {

        double normal1 = GAUSS_RANDOM(1.0),
	  normal2 = GAUSS_RANDOM(1.0),
	  c11 = q[0].cutoff.constant,
	  c12 = q[1].cutoff.constant,
	  c22 = q[3].cutoff.constant,
	  x[2];
	
        if (c11 < 0.0 || c22*c11 - c12*c12 < 0.0 )
	  RFERROR("Cannot simulate field with cutoff, matrix of constants is not positive definite.");
        x[0] = SQRT(c11)*normal1 ;
        x[1] = c12/SQRT(c11)*normal1 + SQRT(c22 - c12*c12/c11 )*normal2;

        if (global->general.vdim_close_together) {
            //one location two values
            for (int i = 0; i < totpts; i++)  {
                for (int j = 0; j < vdim; j++)  {
                    res[i*vdim + j] += x[j];
                }
            }
        } else {
            // the first field, the second field
            for (int j = 0; j < vdim; j++)  {
                for (int i = 0; i < totpts; i++)  {
                    res[i+j*totpts] += x[j];
                }
            }

        }
    }

}

void kappa_localproc(int i, model *cov, int *nr, int *nc){
  kappaGProc(i, cov, nr, nc);
  if (i == CE_MMIN) *nr = 0;

  // $(proj=1) ->intrinsicinternal as $ copies the parameters
  // to intrinsicinternal and at the same time proj=1 reduces
  // the logicaldim
}

void range_co_proc(model VARIABLE_IS_NOT_USED *cov, int k, int i, int j,
		   simple_range_type *range) {
  switch(k) {
  case LOCPROC_DIAM : 
    range->min = 0.0; //  CUTOFF_DIAM
    range->max = RF_INF;
    range->pmin = 1e-10;
    range->pmax = 1e10;
    range->openmin = true;
    range->openmax = true;
    break;
  case LOCPROC_A : 
    range->min = 0.0; // cutoff_a
    range->max = RF_INF;
    range->pmin = 0.5;
    range->pmax = 2.0;
    range->openmin = true;
    range->openmax = true;
    break;
  default : range_ce(cov, k, i, j, range);
  }
}




void range_intrinCE(model VARIABLE_IS_NOT_USED *cov, int k, int i, int j,
		    simple_range_type *range) {
  switch(k) {
  case LOCPROC_DIAM : 
    range->min = 0.0; //  CUTOFF_DIAM
    range->max = RF_INF;
    range->pmin = 0.01;
    range->pmax = 100;
    range->openmin = true;
    range->openmax = true;
    break;
  case LOCPROC_R : 
    range->min = 1.0; // stein_r
    range->max = RF_INF;
    range->pmin = 1.0;
    range->pmax = 20.0;
    range->openmin = false;
    range->openmax = true;
    break;
  default : range_ce(cov, k, i, j, range);
  }
}


void do_circ_embed_intr(model *cov, gen_storage *S) {
  // printf("BASE = %d \n", cov->base == NULL);
  assert(cov->zaehler >= 0);
  model *key = cov->key,
    *sub = key;
  assert(key != NULL);
  sub = sub->key != NULL ? sub->key : sub->sub[0];
     
  location_type *loc = Loc(cov);
  double x[MAXCEDIM], dx[MAXCEDIM],  
    *res = cov->rf;
  Long index[MAXCEDIM];
  int 
    dim = ANYDIM,
    row = dim,
    col = dim,
    rowcol =  row * col;
  GETSTORAGE(s, sub ,   localCE); 
  assert(s != NULL);
  // PMI(cov->calling);
  // printf("%d %d %d %d %d\n", s==NULL, s->correction == NULL, cov->SlocalCE == NULL, cov->calling->SlocalCE == NULL, cov->zaehler);

 #ifdef SCHLATHERS_MACHUINE
  assert(s->correction[0] != 0);
#endif  
  SAVE_GAUSS_TRAFO;

  do_circ_embed(key, S);

  for (int k=0; k<row; k++) {
    index[k] = 0;
    dx[k] = x[k] = 0.0;
  }

  //  PMI(cov->calling);
  //printf("%d\n", cov->key->SlocalCE != NULL);
  if (loc->caniso == NULL) {
    for (int k=0; k<row; k++) {
      double normal = GAUSS_RANDOM(1.0);
      //    printf("k=%d %d %f %f\n", k, row, dx[k], s->correction[0]);
      dx[k] += s->correction[0] * normal; //
    }
  } else {
    for (int i=0; i<rowcol; ) {
      double normal = GAUSS_RANDOM(1.0);
      for (int k=0; k<row; k++)
	dx[k] += s->correction[i++] * normal;
    }
  }
  for (int k=0; k<row; k++) dx[k] *= loc->xgr[k][XSTEP];
  int r;
  for(r=0; ; ) {
    for (int k=0; k<row; k++) res[r] += (double) x[k]; 
    r++;
    int k=0;
    while( (k<row) && (++index[k]>= (int) loc->xgr[k][XLENGTH])) { // OK
      index[k]=0;
      x[k] = 0.0;
      k++;
    }
    if (k>=row) break;
    x[k] += dx[k];
  }
  BOXCOX_INVERSE;
}




int struct_ce_approx(model *cov, model **newmodel) {
  assert(COVNR==CE_CUTOFFPROC_INTERN || COVNR==CE_INTRINPROC_INTERN 
	 || COVNR==CIRCEMBED);
  model *next = cov->sub[0];
  bool cutoff = COVNR ==  CE_CUTOFFPROC_INTERN;
  if (next->pref[COVNR==CIRCEMBED ? CircEmbed : 
		 cutoff ? CircEmbedCutoff : CircEmbedIntrinsic] == PREF_NONE)
    RETURN_ERR(ERRORPREFNONE);
 
  //PMI(cov);


  // assert(!Locgrid(cov));
 
  if (!Locgrid(cov)) {    
    ASSERT_NEWMODEL_NULL;
    location_type *loc = Loc(cov);  
    double max[MAXCEDIM], min[MAXCEDIM],  centre[MAXCEDIM], 
      x[3 * MAXCEDIM],
      approx_gridstep = P0(CE_APPROXSTEP);
    int k, d, len, err,
      maxgridsize = P0INT(CE_APPROXMAXGRID),
      spatialdim = loc->spatialdim;
     
    if (approx_gridstep < 0)
      SERR("approx_step < 0 forbids approximative circulant embedding");
    
    if (OWNTOTALXDIM  != Loctsdim(cov))
      SERR("the dimensions of the coordinates and those of the process differ");
    
    GetDiameter(loc, min, max, centre);
  
    if (loc->Time) {
      if ((int) loc->T[XLENGTH] > maxgridsize) SERR("temporal grid too large");
      maxgridsize /= (int) loc->T[XLENGTH];
    }

    if (ISNAN(approx_gridstep)) {
      // hier assert(false);
      double size = loc->spatialtotalpoints * 3;
      if (size > maxgridsize) size = maxgridsize;
      for (k=d=0; d<spatialdim; d++, k+=3) {
	x[k+XSTART] = min[d];
	x[k+XLENGTH] = (int) POW(size, 1.0 / (double) spatialdim);
	x[k+XSTEP] = (max[d] - min[d]) / (x[k+XLENGTH] - 1.0); // OK
      }
    } else {
     len = 1;
      for (k=d=0; d<spatialdim; d++, k+=3) {
	x[k+XSTART] = min[d];
	len *= 
	  x[k+XLENGTH] = (int) ((max[d] - min[d]) / approx_gridstep)+1;
	x[k+XSTEP] = approx_gridstep;
      }
      if (len > maxgridsize) SERR("userdefined, approximate grid is too large");
    }
  
    if (cov->key!=NULL) COV_DELETE(&(cov->key), cov);
    err = covcpyX(&(cov->key), cov, x, loc->T, spatialdim, spatialdim,
		  UNKNOWN_NUMBER_GRIDPTS, loc->Time, true, false);
    if (err != NOERROR) RETURN_ERR(err);
    SET_CALLING(cov->key, cov);

    if ((err = CHECK_GEN(cov->key,// cov  11.1.19			 
			 VDIM0, VDIM1, GaussMethodType, false))
	//C H E C K(cov, cov->tsdim, cov->xdimprev, cov->typus,
	//	     cov->domprev, cov->isoprev, VDIM0, cov->fr  ame)  
	!= NOERROR) RETURN_ERR(err);
  }

  if (COVNR == CIRCEMBED) { RETURN_NOERROR; }
  else return struct_ce_local(Locgrid(cov) ? cov : cov->key, newmodel);

}




int init_ce_approx(model *cov, gen_storage *S) {  
  if (Locgrid(cov)) {
    if (COVNR==CIRCEMBED) return init_circ_embed(cov, S);
    else return init_circ_embed_local(cov, S); 
  }

  ASSERT_UNREDUCED;

  location_type *loc = Loc(cov),
    *keyloc = Loc(cov->key);  
  assert(cov->key != NULL && keyloc != NULL);
 Long i, cumgridlen[MAXCEDIM], 
   totspatialpts = loc->spatialtotalpoints;
  int d, err,
    // maxgridsize = P0INT(CE_APPROXMAXGRID),
    spatialdim = loc->spatialdim,
   dim = ANYDIM;

  ASSERT_ONESYSTEM;
 
  cov->method = COVNR==CIRCEMBED ? CircEmbed 
    : COVNR== CE_CUTOFFPROC_INTERN ? CircEmbedCutoff : CircEmbedIntrinsic;

  if (loc->distances) RETURN_ERR(ERRORFAILED);
  NEW_STORAGE(approxCE);
  ALLC_NEWINT(SapproxCE, idx, totspatialpts, idx);
  
  cumgridlen[0] = 1;
  for (d=1; d<dim; d++) {
    cumgridlen[d] =  cumgridlen[d-1] * (int) keyloc->xgr[d-1][XLENGTH];    
  }

  double *xx = loc->x;

  for (i=0; i<totspatialpts; i++) {
    int dummy;    
    for (dummy = d = 0; d<spatialdim; d++, xx++) {
      dummy += cumgridlen[d] *
	(int) ROUND((*xx - keyloc->xgr[d][XSTART]) / keyloc->xgr[d][XSTEP]);
    }    
    idx[i] = dummy;
    assert(dummy >= 0);
  }

  if (COVNR==CIRCEMBED) err = init_circ_embed(cov->key, S);
  else err = init_circ_embed_local(cov->key, S);
  if (err != NOERROR) RETURN_ERR(err);
  
  if ((err = ReturnOwnField(cov)) != NOERROR) RETURN_ERR(err);
 
  cov->key->initialised = cov->simu.active = err == NOERROR;
  RETURN_NOERROR;
}

void do_ce_approx(model *cov, gen_storage *S){
  if (Locgrid(cov)) {
    if (COVNR==CIRCEMBED) do_circ_embed(cov, S);
    else if (COVNR== CE_CUTOFFPROC_INTERN) do_circ_embed_cutoff(cov, S);
    else do_circ_embed_intr(cov, S);
    return;
  }

  
  //  PMI(cov);

  model *key=cov->key;
  location_type *loc = Loc(cov);
  getStorage(s ,   approxCE); 
  //model *cov = meth->cov;
  Long i;
  int
     vdim = VDIM0,
    *idx = s->idx;
  double 
    *res = cov->rf,
     *internalres = key->rf;

  // APMI(key);

  DO(key, S);
  location_type *key_loc = Loc(key);
  if (key_loc->Time) { // time separately given   
    Long  t,
      j = 0,
      instances = (int) loc->T[XLENGTH],
      totspatialpts = loc->spatialtotalpoints,
      gridpoints = Loc(key)->spatialtotalpoints;
    
    for (int v=0; v<vdim; v++){
      for (t=0; t<instances; t++, internalres += gridpoints) {
	for (i=0; i<totspatialpts; i++) {
	  res[j++] = internalres[idx[i]];
	}
      }
    }
  } else { // no time separately given
    Long totpts = loc->totalpoints;
    int 
      j = 0,
      totalkey = Loc(key)->totalpoints;
    for (int v=0; v<vdim; v++, internalres += totalkey){
      for (i=0; i<totpts; i++) {
      //    printf("i=%d %10g %d\n", i, loc->x[i], idx[i]);
      // print(" %d %d %10g %10g", i, idx[i], loc->x[i*2], loc->x[i*2 +1]);
      // print("    %10g\n", internalres[idx[i]]);
	res[j++] = internalres[idx[i]];
      }
    }
  }

}

