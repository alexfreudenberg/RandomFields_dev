/*
Authors
Martin Schlather, schlather@math.uni-mannheim.de

main library for unconditional simulation of random fields

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
  
 
#include <Rmath.h>  
#include <stdio.h>  
#include <string.h> 
#include "questions.h"
#include "variogramAndCo.h"

model* wheregenuineStatOwn(model *cov) {
  model *sub = cov;
  if (equalsnowGaussMethod(sub) || SUBNR==GAUSSPROC) {
    sub = sub->sub[0];
    while (equalsnowGaussMethod(sub) || SUBNR==GAUSSPROC) sub = sub->sub[0];
  } else if (isnowProcess(sub)) {
    NotProgrammedYet("");
  };

  if (cov->pref[Nothing] == PREF_NONE ||
      !(isnowPosDef(sub) || isXonly(SYSOF(sub)))) {
    // Variogramme sind hier definitiv erlaubt, da durch Addition einer 
    // Konstanten, das Zeug zu einer Kovarianzmatrix gemacht werden kann
    // siehe direct.cc    
    //assert(({PMI(cov, "cov matrix"); true;})); //
    //PMI(cov);
    ERR("covariance matrix: given model is not a covariance function");
  }
   
  return sub;
} 


#define STANDARDSTART(ignore_x, ignore_y, TOTX, TOTY, I_BASE)	\
  model *cov = Cov;					\
  assert(cov != NULL); 				\
  if (equalsnowGaussMethod(cov) || COVNR==GAUSSPROC) cov = cov->sub[0];	\
  model *calling = cov;				\
  if (calling->Sfctn == NULL) {   \
     assert(isnowVariogram(calling)); 		\
    calling = cov->calling; /* either interface or process	*/	\
    if (calling != NULL && calling->Sfctn == NULL) {			\
      assert(isnowProcess(calling));					\
      calling = calling->calling;					\
    }									\
  }									\
  assert(calling != NULL);						\
  if ((calling)->Sfctn == NULL) { BUG; }				\
  model *genuine = wheregenuineStatOwn(cov); /* might have gaussmethod-proc 28.12.20: genuine nicht durch cov ersetzten!! (Konflikt in FINISH_START) */;\
  FINISH_START(calling, cov, ignore_x, ignore_y, 0, TOTX, TOTY, I_BASE); \
  double *cross = fctn->cross			        

//  printf("OKcc\n\n");PMI0(calling);PMI0(cov);assert(calling != NULL && (equalsnowInterface(calling) || isnowProcess(calling))); printf("OK\n\n");


/* assert(VDIM0 == VDIM1); */	
void CovVario(model *Cov, bool is_cov, bool pseudo, int select, bool ignore_y,
	      double *v);

void Variogram(model *cov, double *v) {
  CovVario(cov, false, false, NA_INTEGER, false, v);
}
void Covariance(model *cov, double *v) {
  CovVario(cov, true, false, NA_INTEGER, false, v);
}
void CovarianceT(model *cov, int base_i_row, double *v) {
  CovVario(cov, true, false, base_i_row, true, v);
  // ignore_y=true : swap of x and y!
}
void Pseudomadogram(model *cov, double alpha, double *v) {
  model *calling = cov->calling;
  CovVario(cov, false, true, NA_INTEGER, false, v);
  if (alpha < 2.0) {
    Long tot = Loctotalpoints(cov),
      vdim = VDIM0,      
      end = vdim * vdim * tot;
    double half = 0.5 * alpha,
       c = POW(2, alpha - 1) * gammafn(0.5 + half) / SQRTPI;
    if (alpha == 1.0) for (int i=0; i < end; i++) v[i] = c * SQRT(v[i]);
    else for (int i=0; i < end; i++) v[i] = c * POW(v[i], half);
  }
}

void CovVario(model *Cov, bool is_cov, bool pseudo, int select, bool ignore_y,
	      double *v) {
  globalparam *global = &(Cov->base->global);
  // note: here, only a vector is return (in the univatiate case),
  //       not a matrix. So, distances do not make sense.
  // if kernel and y not given then y:=0
  // if x and y are given, length of x is taken und y is recycled
 
  STANDARDSTART(ignore_y, ignore_y, 0, 1, // 12.1.21 select!=NA_INTEGER,
		select==NA_INTEGER ? 0 : select);
  info[INFO_EXTRA_DATA_X] = ignore_y;
  info[INFO_EXTRA_DATA_Y] = !ignore_y;
     
   double 
    *zero = ZERO(cov),
   *C0x = fctn->C0x,
    *C0y = fctn->C0y;
  int vdimP1 = vdim0 + 1;

  assert(cov != NULL);
  assert(cov->base != NULL);
  assert(cov->ownloc != NULL || cov->prevloc != NULL);
  //printf("set = %d %d %d %ld %ld\n", cov->base->set,cov->ownloc != NULL, cov->prevloc != NULL, cov->ownloc,cov->prevloc  );
  //PMI0(cov);
  assert(LocLoc((cov)->prevloc, cov) != NULL);
  
  assert(Loc(cov) != NULL);
  if (LocDist(cov)) BUG;

  bool kernel = equalsKernel(DOM(PREVSYSOF(genuine), 0));

  if (kernel && !ygiven && PL > 0 &&
      cov->base->global.messages.warn_singlevariab){
    WARN1("'%.50s' is called with a single variable only, although it is used as a kernel. So, the second variable is set to zero, here.\n", NICK(cov));
  }
 

  if (is_cov) {
    assert(({/*P M I(cov, "Cov/vario"); */ cov->pref[Nothing] != PREF_NONE &&
	    isnowShape(genuine);})); //
  } else {
    bool isvario =  isnowVariogram(genuine);
   //printf("hier !!\n");
    if (cov->pref[Nothing] == PREF_NONE || !isvario) {
      assert(({PMI(cov); true;})); //
      ERR("given model is not a variogram");
    }
    if (!kernel && isCartesian(PREV)) {   
      if (!isCartesian(OWN)) BUG;
      COV(zero, info, genuine, C0y);
      //printf("C0y=%f\n", *C0y);
    } else {
      if (vdim0 > 1 && kernel) 
	ERR("multivariate variogram only calculable for stationary models");
      NONSTATCOV(zero, zero, info, genuine, C0y);
    }
    if (pseudo && vdim0 > 1) { 
      for (int m=0; m<vdim0; m++) 
	for (int n=0; n<=m; n++) 
	  C0y[n * vdim0 + m] = C0y[m * vdim0 + n] =
	    0.5 * (C0y[n * vdimP1] + C0y[m * vdimP1]); // mean of the variances
    }
  }
    

#define UNIVAR COV(x, info, genuine, v)
#define UNIVAR_Y NONSTATCOV(x, y, info, genuine, v) 

#define VARIO_UNIVAR	  \
  COV(x, info, genuine, cross);	  \
  *v = *C0y - *cross; 

#define VARIO_UNIVAR_Y	   	\
  NONSTATCOV(x, y, info, genuine, cross);		\
  NONSTATCOV(x, x, info, genuine, C0x);		\
  NONSTATCOV(y, y, info, genuine, C0y);		\
    *v = 0.5 * (*C0x + *C0y) - *cross;

#define MULT 					\
  COV(x, info, genuine, cross);				\
  VDIM_LOOP(cross[u]);

#define MULT_Y		\
  NONSTATCOV(x, y, info, genuine, cross);	\
  VDIM_LOOP(cross[u])
  
#define PSEUDO_MULT 							\
  COV(x, info, genuine, cross);						\
  VDIM_LOOP(C0y[u] - cross[u])

#define VARIO_MULT 							\
  COV(x, info, genuine, cross);						\
  VDIM_LOOP(C0y[u] - 0.5 * (cross[u] + cross[w]))	
 
#define VARIO_MULT_Y 							\
  NONSTATCOV(x, y, info, genuine, cross);				\
  NONSTATCOV(x, x, info, genuine, C0x);					\
  NONSTATCOV(y, y, info, genuine, C0y);					\
  VDIM_LOOP(0.5 * (C0x[u] + C0y[u] - cross[u] - cross[w]))
 
#define PSEUDO_MULT_Y							\
  NONSTATCOV(x, y, info, genuine, cross); 				\
  NONSTATCOV(x, x, info, genuine, C0x);					\
  NONSTATCOV(y, y, info, genuine, C0y);					\
  double *C = v + VDIM_0 * i_row + i_col * VDIMtotX;			\
  int u = 0;								\
  int m = 0;								\
  for (Long n1=0; m<vdimSq; m+=vdimP1, n1+=NINCR) {			\
    int w = 0;								\
    for (Long m1=n1; w<vdimSq; m1+=MINCR, w+=vdimP1){		\
      C[m1] = 0.5 * (C0x[m] + C0y[w]) - cross[u++];			\
    }									\
  }
      
  // printf("covvario 2 %.50s, line %d %d %ld %ld %ld\n",				
  //     __FILE__, __LINE__, ygiven, y, zero, fctn->y) ;	

  //  printf("grid = %d ygiven=%d kernel=%d trafo=%d\n", grid, ygiven, kernel, trafo);
  
  assert(y == fctn->y);
  if (select != NA_INTEGER) {
    if (caniso != NULL) BUG;
    if (trafo) // TO DO --- should not appear as ExpandGrid in kriging.R
      ERR("kriging and conditional simulation for spatio-temporal fields are not programmed yet");
    // trafoY wird hier komplett umgangen:
    if (gridY) { // works for both grid and !grid (for x)
      int r = i_col_base;
 	assert(y == fctn->y);
     for (d=0; d<tsxdim; d++) {
	int n = r % (int) grY[d][XLENGTH];
	r /= grY[d][XLENGTH];	
	incy[d] = grY[d][XSTEP];	   
	y[d] = ystart[d] = grY[d][XSTART] + incy[d] * n;
	endy[d] = 1;
	ny[d] = startny[d] = 0; 
      }
    } else {
      int spatialdim = Locspatialdim(cov),
	spptsY = LocspatialpointsY(cov);
      MEMCOPY(y, LocY(cov) + spatialdim * (i_col_base % spptsY),
	      spatialdim * sizeof(double));
      if (Time) {
	double *T = LocTY(cov);
	y[spatialdim] = T[XSTART] + (i_col_base / spptsY) * T[XSTEP];
      }
      for (d=0; d<tsxdim; d++) { // works for grid (x)
	incy[d] = 1.0;	   
	ystart[d] = y[d];
	endy[d] = 1;
	ny[d] = startny[d] = 0; 
      }
    }
  } else {
    PERFORM_PREPARE;
  }

  //  printf("covvario ygien %d %d tot=%ld %ld\n", ygiven, kernel, totX, totY);
  
  //   printf("grid = %d %d %d\n", grid, ygiven, kernel);
  //    PMI0(cov);
  //if (ygiven && !kernel) crash();

    if (is_cov) {  
    PERFORM(UNIVAR, MULT, UNIVAR_Y, MULT_Y);
  } else if (pseudo) {    
    PERFORM(VARIO_UNIVAR, PSEUDO_MULT, VARIO_UNIVAR_Y, PSEUDO_MULT_Y);
  } else {
     PERFORM(VARIO_UNIVAR, VARIO_MULT, VARIO_UNIVAR_Y, VARIO_MULT_Y);
  }
  STANDARD_ENDE;
} 
 
 
#define swap(x, y) { double swapdummy = x;  x = y;  y = swapdummy; }
#define swapInt(x, y) { int swapdummy = x;  x = y;  y = swapdummy; }

void CovarianceMatrix(model *Cov, bool ignore_y, double *v) {
  
  // NOTE: if LocHasY & !ignore_y then LocY is taken and not LocX

#define  StartCovarianceMatrix(TOTX, TOTY)			\
  globalparam *global = &(Cov->base->global);			\
  STANDARDSTART(!ignore_y, ignore_y, TOTX, TOTY, 0);		\
  bool dist = LocDist(cov);  /* LocxdimOZ(cov) == 1 nicht notwendig! */ \
  Long totM1 = totX - 1;					\
  double *x0 = NULL;	/* never free it  */			\
  info[INFO_EXTRA_DATA_X] = info[INFO_EXTRA_DATA_Y] = !ignore_y
  

  StartCovarianceMatrix(0, 0);
  
#define MULTICOV						\
  double *C = v + VDIM_0 * i_row + i_col * VDIMtotX,		\
    *D = v + VDIM_0 * i_col + i_row * VDIMtotX;			\
  int l=0;							\
  for (Long n=0, m0=0; n<NEND; n+=NINCR, m0+=MINCR) {		\
    Long endfor = n + ENDFORINCR;				\
    for (Long m=n, n0=m0; m<endfor; m+=MINCR, n0+=NINCR) {	\
      C[m] = D[n0] =cross[l++];					\
    }								\
  }								       
  
  bool kernel = equalsKernel(DOM(PREVSYSOF(genuine), 0));
  assert(y == fctn->y);

  if (grid) {
    for (d=0; d<tsxdim; d++){		
      incy[d] = gr[d][XSTEP];
      endy[d] = gr[d][XLENGTH];
      startny[d] = 0;	
    }
    if (!kernel) { // i.e. stationary, i.e. in the case of tsdim=1 
      // we have the same value on any (off)diagonal. Similar for tsdim>1.
      // So using this fact, algorithm gets much faster.
      int lastD = tsxdim - 1,
	*cum = (int*) fctn->cum,
	k = 0;
      cum[0] = 1;
      for (d=1; d<=lastD; d++) cum[d] = cum[d-1] * end[d-1];//==endy
      while (true) {
	i_col = i_row;
	info[INFO_IDX_Y] = i_col;
	for (d=0; d<tsxdim; d++) {
	  y[d] = x[d];
	  ystart[d] = xstart[d];
	  ny[d] = nx[d];
	}
	while (true) {
	  if ( (nx[k] == 0L || ny[k] == 0L))// previous k good guess for current
	    for (k=0 ; k<=lastD; k++) if (nx[k] > 0L &&  ny[k] > 0L) break;
	  if (k > lastD) {
	    NONSTAT2STATCOV(x, y, info, genuine, cross);
	    MULTICOV;
	  } else {
	    int i_col_alt = i_col - cum[k],
	      i_row_alt = i_row - cum[k];
	    assert(i_col_alt >=0 && i_row_alt >= 0);
	    double *Calt = v + VDIM_0 * i_col_alt + i_row_alt * VDIMtotX,  
	      *C = v + VDIM_0 * i_col + i_row * VDIMtotX;

	    for (int n=0; n<NEND; n+=NINCR) {
	      Long endfor = n + ENDFORINCR;	
	      for (int m=n; m<endfor; m+=MINCR) {				
		C[m] = Calt[m];						
	      }								
	    }								
	    								
	    if (i_col != i_row) {					
	      C = v + VDIM_0 * i_row + i_col * VDIMtotX;		
	      Calt = v + VDIM_0 * i_row_alt + i_col_alt * VDIMtotX;	
	      for (int m=0; m<ENDFORINCR; m+=MINCR) {			
		for (int n=m; n<NEND; n+=NINCR) {				
		  C[n] = Calt[n];					
		}							
	      }								
	    }
	    
	  }
	  STANDARDINKREMENT_Y;
	  //	APMI(genuine);
	  if (d >= tsxdim) break; 
	}
	STANDARDINKREMENT_X;
      }
      
    } else { // KERNEL
      while (true) {
	i_col = i_row;
	info[INFO_IDX_Y] = i_col;
	
	for (d=0; d<tsxdim; d++) {
	  y[d] = x[d];
	  ystart[d] = xstart[d];
	  ny[d] = nx[d];
	}
	while (true) {
	  NONSTATCOV(x, y, info, genuine, cross);
	  MULTICOV; 
	  STANDARDINKREMENT_Y;
	  //	APMI(genuine);
	  if (d >= tsxdim) break; 
	}
	STANDARDINKREMENT_X;
      }
    }
    
  } else { // not a grid
    if (trafo) {
      localdim = TransformLoc(cov, NULL, ygiven ? NULL : &xx,
			      ygiven ? &xx : NULL,
			      False);
      assert(localdim == tsxdim);
      x0 = xx;
    } else x0 = LocY(cov, ignore_y);
    x = x0;

    // i_row/col fkt nicht mit parallel!!!!
    // am saubersten waere es, das i als zusaetzliche Koordinate
    // in einem extra koordinatensystem zu uebertragen
    // auch muss nugget zu non-stat nugget umgeschrieben werden,
    // so dass nur fuer nicht-stat. cov-fkten die weitere Koordinate
    // uebertragen werden muss + zusaetzlich die Abfrage "stationaer"
    // ausreicht, um loc->i* auszuschliessen. Oder aber ein Flag
    // setzen in model, das anzeigt, ob loc->i in einem untermodel
    // vorliegt Auch nugget (messfehler) ist eigentlich nicht-stationaer!
    if (dist) {
      double zero = 0.0;      
      for (i_row=0; i_row<totX; i_row++) {
	info[INFO_IDX_X] = i_row;
	x = x0 + localdim * i_row;
	for (y=x, i_col=i_row; i_col<totX; i_col++, y+=localdim) {
	  info[INFO_IDX_Y] = i_col;
	  if (i_col==i_row) COV(&zero, info, genuine, cross)
	    else COV(x0 + (i_row * totM1 - (i_row * (i_row + 1)) / 2 + i_col -1)
		     * localdim, info, genuine, cross);
	  MULTICOV; 
	}
      }
    } else { // no distances
      for (i_row=0; i_row<totX; i_row++) {
	x = x0 + localdim * i_row;
	info[INFO_IDX_X] = i_row;
	for (y=x, i_col=i_row; i_col<totX; i_col++, y+=localdim) {
	  //	  printf("%f %f; %f %f  %d\n", x[0], x[1], y[0], y[1], kernel);
	  info[INFO_IDX_Y] = i_col;

	  //PMI(genuine);	  printf("kernel =  %d\n", kernel);
	  
	  if (kernel) NONSTATCOV(x, y, info, genuine, cross)
	  else NONSTAT2STATCOV(x, y, info, genuine, cross);
	  if (!R_FINITE(cross[0])) GERR("model creates non-finte values.");	
	  MULTICOV; 
	}
      }
    } // no distances
  } // not a grid
  

  if (false) {
    for (int m=0; m<27; m++) {
      for (int n=0; n<27; n++) {
	//printf("%+2.2f ", v[n * 27 + m]);
      }
      //printf("\n"); 
    }
  }

 ErrorHandling: 
  STANDARD_ENDE;
  if (err!=NOERROR) XERR(err); 
   
} // CovarianzMatrix


void CovarianceMatrix(model *Cov, bool ignore_y, int *idx, int Nidx,
		      double *v) {
  // only a Nidx x Nidx matrix is calculated with indices given by idx.
  // NOTE: if LocHasY & !ignore_y then LocY is taken and not LocX 
  StartCovarianceMatrix(Nidx, Nidx);

  bool kernel = equalsKernel(DOM(PREVSYSOF(genuine), 0));
  assert(y == fctn->y);
  
  if (grid) {
    for (i_row=0; i_row<Nidx; i_row++) {
      int I = info[INFO_IDX_X] = idx[i_row];
      for (d=0; d<tsxdim; d++) {
	x[d] = gr[d][XSTART] + (I % (int) gr[d][XLENGTH]) * gr[d][XSTEP];
	I /= gr[d][XLENGTH];
      }
      for (i_col = i_row; i_col < Nidx; i_col++) {
	I = info[INFO_IDX_Y] = idx[i_col];
	assert(y == fctn->y);
	for (d=0; d<tsxdim; d++) {
	  y[d] = gr[d][XSTART] + (I % (int) gr[d][XLENGTH]) * gr[d][XSTEP];
	  I /= gr[d][XLENGTH];
	}
	if (kernel) NONSTATCOV(x, y, info, genuine, cross)
	else NONSTAT2STATCOV(x, y, info, genuine, cross);
	if (!R_FINITE(cross[0])) GERR("model creates non-finte values.");	
	MULTICOV;
      }
    }  
  } else { // not a grid
    if (trafo) {
      localdim = TransformLoc(cov, NULL, ygiven ? NULL : &xx,
			      ygiven ? &xx : NULL, False);
      assert(localdim == tsxdim);
      x0 = xx;
    } else x0 = LocY(cov, ignore_y);
    x = x0;

    if (dist) {
      double zero = 0.0;      
      for (i_row=0; i_row<Nidx; i_row++) {
	int I = info[INFO_IDX_X] = idx[i_row];	
	x = x0 + localdim * I;
	for (i_col = i_row; i_col < Nidx; i_col++) {	  
	  int J = info[INFO_IDX_Y] = idx[i_col];
	  y = x0 + localdim * J;
	  if (I==J) COV(&zero, info, genuine, cross)
	  else COV(x0 + (I * totM1 - (I * (I + 1)) / 2 + J -1)
		   * localdim, info, genuine, cross);
	  MULTICOV; 
	}
      }
    } else { // no distances
      int I;
      for (i_row=0; i_row<Nidx; i_row++) {
	I = info[INFO_IDX_X] = idx[i_row];
	x = x0 + localdim * I;
	for (i_col = i_row; i_col < Nidx; i_col++) {
	  int J = info[INFO_IDX_Y] = idx[i_col];
	  y = x0 + localdim * J;
	  if (kernel) NONSTATCOV(x, y, info, genuine, cross)
	  else NONSTAT2STATCOV(x, y, info, genuine, cross);
	  if (!R_FINITE(cross[0])) GERR("model creates non-finte values.");	
	  MULTICOV;
	}
      }
    } // no distances
  } // not a grid
  

  if (false) {
    for (int m=0; m<27; m++) {
      for (int n=0; n<27; n++) {
	//printf("%+2.2f ", v[n * 27 + m]);
      }
      //printf("\n"); 
    }
  }

 ErrorHandling: 
  STANDARD_ENDE;
  if (err!=NOERROR) XERR(err); 
   
} // CovarianzMatrix





void CovarianceMatrixCols(model *Cov, bool ignore_y, int row, double *v) {

  // NOTE: if LocHasY and !ignore_y, then LocY is taken not LocX
  StartCovarianceMatrix(1, 0);
  assert(y == fctn->y);

    
#define MULTICOV_COL							\
  double *C = v + VDIM_0 * 0 + i_col * VDIMtotX;			\
  int l = 0;								\
  for (Long m=0; m<ENDFORINCR; m+=MINCR) {				\
    for (Long n=m; n<NEND; n+=NINCR) {					\
      C[n] = cross[l++];						\
    }									\
  }

  info[INFO_IDX_X] = i_row = row;
  info[INFO_IDX_Y] =i_col = 0;
  if (grid) {
    int r = row;
    for (d=0; d<tsxdim; d++) {
      incy[d] = gr[d][XSTEP];   
      y[d] = ystart[d] = gr[d][XSTART];
      endy[d] = gr[d][XLENGTH]; 
      ny[d] = startny[d]= 0;

      nx[d] = r % end[d];
      r /= end[d];	
      x[d] = ystart[d] + incy[d] * (double) nx[d];
    }
    while (true) {
      NONSTATCOV(x, y, info, genuine, cross);
      MULTICOV_COL; 
      STANDARDINKREMENT_Y;
      //	APMI(genuine);
      if (d >= tsxdim) break; 
    }
  } else { 
    if (trafo) {
      localdim = TransformLoc(cov, NULL, ygiven ? NULL : &xx,
			      ygiven ? &xx : NULL, False);    
      assert(localdim == tsxdim);
      x = xx;
    } else x = LocY(cov, ignore_y);
 
    if (dist) {
      double zero = {0.0};
      for (y=x, i_col=0; i_col<totX;  i_col++, y+=localdim) {
	info[INFO_IDX_Y] = i_col;
	if (i_col == row) {
	  COV(&zero, info, genuine, cross); 
	} else {
	  int i_min = MIN(i_col, row),
	    i_max = MAX(i_col, row);
	  COV(x + (i_min * totM1 - (i_min * (i_min + 1)) / 2 + i_max -1) *
	      localdim, info, genuine, cross);
	}
	MULTICOV_COL;
      }
    } else {// no distances
      x0 = x + localdim * row;      
      for (y=x, i_col=0; i_col<totX; i_col++, y+=localdim) {
	info[INFO_IDX_Y] = i_col;
	NONSTATCOV(x0, y, info, genuine, cross);
	assert(R_FINITE(cross[0]));
  	MULTICOV_COL;
      }
    } // no distances
  } // not a grid
  

  if (false) {
    for (int m=0; m<27; m++) {
      for (int n=0; n<27; n++) {
	//printf("%+2.2f ", v[n * 27 + m]);
      }
      //printf("\n"); 
    }
  }
  //assert(false);
  
  STANDARD_ENDE;
  if (err!=NOERROR) XERR(err); 
  
  //  int i,j,k;
  //  for (k=i=0; i<totX*totX; i+=totX) {
  //    for (j=0; j<totX;j++) printf("%10g ", v[k++]);
  //    printf("\n");  }
  
} // CovarianzMatrixCol


void InverseCovMatrix(model *cov, double *v, double *det) {// currently unused
  // needed when Markov models are implemented
  Long vdimtot = (Long) Loctotalpoints(cov) * VDIM0;
  assert(VDIM0 == VDIM1);
  Covariance(cov, v);
  if (cov->Ssolve == NULL) SOLVE_STORAGE;
  int Exterr = Ext_solvePosDef(v, vdimtot, true, NULL, 0, det, cov->Ssolve);
  if (Exterr != NOERROR){
    STRCPY(cov->err_msg, cov->Ssolve->err_msg);
    OnErrorStop(Exterr, cov->err_msg);
  }
}


//////////////////////////////////////////////////////////////////////
// Schnittstellen

/*
#define STANDARDINTERN					\
  if (reg < 0 || reg > MODEL_MAX) XERR(ERRORREGISTER);	\
   assert(currentNrCov != UNSET);		\
   model *cov = KEY()[reg];				\
  if (cov == NULL) { ERR("register not initialised") }	\
  model *truecov = !equalsnowInterface(cov) ?		\
    cov : cov->key == NULL ? cov->sub[0] : cov->key
*/
 
#define STANDARDINTERN_SEXP_BASIC					\
  if (INTEGER(reg)[0] < 0 || INTEGER(reg)[0] > MODEL_MAX) XERR(ERRORREGISTER); \
   assert(currentNrCov != UNSET);				\
   model *cov = KEY()[INTEGER(reg)[0]];				\
  if (cov == NULL) { ERR("register not initialised") }			

#define STANDARDINTERN_SEXP						\
  STANDARDINTERN_SEXP_BASIC;						\
  model VARIABLE_IS_NOT_USED *truecov =  !equalsnowInterface(cov)	? \
    cov : cov->key == NULL ? cov->sub[0] : cov->key;			\
  if (equalsnowGaussMethod(truecov)) truecov = truecov->sub[0]

//  if (cov->pref[Nothing] == PREF_NONE) { PMI(cov); XERR(ERRORINVALIDMODEL) }


SEXP CovLocNonGrid(SEXP reg, SEXP x, SEXP Y, SEXP result) { 
   STANDARDINTERN_SEXP;
  double *y = TYPEOF(Y)==NILSXP ? NULL : REAL(Y);
  int err,
    dim = nrows(x),
    lxy = ncols(x);
  raw_type save = cov->base->rawConcerns;
 
  assert(dim == Loctsdim(cov));
  assert(y == NULL || (ncols(Y) == lxy && nrows(Y) == dim));

  if (LocSets(cov) > 1) BUG;

  
  location_type **old = cov->prevloc;
  assert(cov->ownloc == NULL);
  LOC_DELETE(&(cov->prevloc));
  if ((err = partial_loc_set(Loc(cov), REAL(x), y, lxy, y == NULL ? 0 : lxy,
			     false, dim, NULL, NULL, false, false, false))
      != NOERROR) XERR(err);   
  SetLoc2NewLoc(cov, old, cov->prevloc,
#ifdef SCHLATHERS_MACHINE
	   true
#else
	   false
#endif
	   );

  cov->base->rawConcerns = neverRaw;
  Covariance(truecov, REAL(result));
  cov->base->rawConcerns = save;

  location_type *loc = Loc(cov);
  loc->totalpoints = loc->totalpointsY = 0; // ok
  loc->x = loc->Y = NULL; // ok

  return R_NilValue;
}


SEXP LocNonGrid(SEXP reg, SEXP x) { 
   STANDARDINTERN_SEXP;
  int err;
#ifdef SCHLATHERS_MACHINE
  int tsdim = Loctsdim(cov),
    xdimOZ = LocxdimOZ(cov),
    spatialdim = Locspatialdim(cov);
  bool Time = LocTime(cov);
#endif

  location_type **old = cov->prevloc,
    **neu = loc_set(x);
  SetLoc2NewLoc(cov, old, neu,
#ifdef SCHLATHERS_MACHINE
	   true
#else
	   false
#endif
	   );
  assert(cov->prevloc == neu);
  LOC_DELETE(&old); 
  cov->base->rawConcerns = neverRaw;
 
#ifdef SCHLATHERS_MACHINE
  assert(tsdim == Loctsdim(cov));
  assert(xdimOZ == LocxdimOZ(cov));
  assert(spatialdim == Locspatialdim(cov));
  assert(Time == LocTime(cov));
#endif  

  return R_NilValue;
}



SEXP MomentsIntern(SEXP reg, SEXP Alpha) { // used
  STANDARDINTERN_SEXP;
  SEXP ans;
  int 
    vdim =VDIM0,
    size = Loctotalpoints(cov) * vdim * vdim;
  double alpha = REAL(Alpha)[0];
  PROTECT(ans = allocVector(REALSXP, size));
  if (alpha == VARIOGRAM) Variogram(truecov, REAL(ans));
  else if (alpha == COVARIANCE) Covariance(truecov, REAL(ans));
  else Pseudomadogram(truecov, FABS(alpha) /* safe */, REAL(ans)); 
  UNPROTECT(1);
  return(ans);
}

