

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017  Martin Schlather

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



#ifndef variogramAndCo_H
#define variogramAndCo_H 1

#define STANDARDSTART_X							\
  assert(loc != NULL && loc->xgr[0] != NULL);				\
  for (d=0; d<tsxdim; d++) {						\
    inc[d] = loc->xgr[d][XSTEP];					\
    start[d] = 0;							\
    end[d] = loc->xgr[d][XLENGTH];					\
    nx[d] = start[d];							\
    x[d] = xstart[d] = loc->xgr[d][XSTART];				\
  }									\
  x[tsxdim] = 0



#define FINISH_START_X(whereSfctn)					\
  fctn_storage *fctn= (whereSfctn)->Sfctn;				\
  if (fctn == NULL) BUG;						\
  location_type *loc = Loc(whereSfctn);					\
  assert(loc != NULL);							\
  /*bool grid =loc->grid && !loc->Time && loc->caniso==NULL, why !time ?? 10.3.15*/		\
  bool grid =loc->grid && loc->caniso==NULL,				\
    trafo = (loc->Time || loc->caniso != NULL) && !loc->distances;	\
  Long tot = loc->totalpoints;						\
  int d,								\
    *start=fctn->start,							\
    *end=fctn->end,							\
    i_row = 0,								\
    *nx=fctn->nx,							\
    tsxdim = PREVTOTALXDIM; /* weicht u.U. von logicaldim bei dist=true ab */ \
  if (end == NULL) BUG;							\
  double *x = fctn->x,							\
    *xstart= fctn->xstart,						\
    *inc=fctn->inc,							\
    *xx = NULL;		/* dito			   */			\
  assert(x != NULL);							\
  if (grid) {								\
    STANDARDSTART_X;							\
  }									



#define FINISH_START(whereSfctn, whereVdim, ignore_y, vdim2)		\
  FINISH_START_X(whereSfctn);		\
  /*bool grid =loc->grid && !loc->Time && loc->caniso==NULL, why !time ?? 10.3.15*/		\
  bool ygiven = !(ignore_y) && (loc->y != NULL || loc->ygr[0] != NULL);/* might be changed !*/ \
  int i_col = 0,							\
    vdim0  = (whereVdim)->vdim[0],					\
    vdimSq = vdim0 * (whereVdim)->vdim[vdim2], /*not nec. squ!!*/	\
    *endy =fctn->endy,							\
    *startny =fctn->startny,						\
    *ny = fctn->ny;	 						\
 Long VARIABLE_IS_NOT_USED vdim0tot = vdim0 * tot,			\
    vdimSqtot = vdimSq * tot,						\
    err = NOERROR;							\
  double *ystart =fctn->ystart,						\
     *yy = NULL,		/* set by Transform Loc; freed  */	\
     *incy =fctn->incy			
  


#define INCLUDE_VAL							\
  double **Val= fctn->Val,						\
    *cross=fctn->cross;							\
  if (vdimSq > 1) {							\
    for (Long v=0, i = 0; i<vdimSqtot ; i+=vdim0tot) {			\
      for (Long ii=i, j=0; j<vdim0; ii+=tot, v++, j++) {		\
	Val[v] = value + ii;						\
      }									\
    }									\
  }									\
  if (!(kernel && !ygiven && PL > 0 && GLOBAL.internal.warn_singlevariab)) { } else {			\
      WARN1("'%.50s' is called with a single variable only, although it is used as a kernel. So, the second variable is set to zero, here.\n", NICK(cov));  \
  } 



#define STANDARDINKREMENT_X	\
  d = 0;			\
  nx[d]++;			\
  x[d] += inc[d];		\
  while (nx[d] >= end[d]) {	\
    nx[d] = start[d];		\
    x[d] = xstart[d];		\
    if (++d >= tsxdim) break;	\
    nx[d]++;			\
    x[d] += inc[d];		\
  }				\
  x[tsxdim] = ++i_row;		\
  if (d < tsxdim) { } else break		
  // end define StandardInkrement


				
// Achtung hier *nicht* incy, starty
#define STANDARDINKREMENT_Y			\
  d = 0;					\
  ny[d]++;					\
  y[d] += incy[d];				\
  while (ny[d] >= endy[d]) {			\
     ny[d] = startny[d];			\
     y[d] = ystart[d];				\
     if (++d >= tsxdim) break;			\
     ny[d]++;					\
    y[d] += incy[d];				\
   }						\
   y[tsxdim] = ++i_col				\
   // end define Stand ardInkrement_Y
                 

// ny wird nur inv CovMatrix verwendet! Stattdessen dort weder
// ystart noch incy!
#define STANDARDSTART_Y_SUPPL						\
  for (d=0; d<tsxdim; d++){						\
    incy[d] = loc->ygr[d]==NULL ? loc->xgr[d][XSTEP] : loc->ygr[d][XSTEP]; \
    y[d] = ystart[d] =							\
      loc->ygr[d]==NULL ? loc->xgr[d][XSTART] : loc->ygr[d][XSTART];	\
    endy[d] = loc->ygr[d]==NULL ? loc->xgr[d][XLENGTH]:loc->ygr[d][XLENGTH]; \
    startny[d] = 0;							\
    ny[d] = startny[d];							\
  }									\
  y[tsxdim] = i_col = 0
								

#define RECYCLE_Y				\
  if (d <= tsxdim) { } else {			\
    for (d=0; d<tsxdim; d++) {			\
      y[d] = ystart[d] = loc->ygr[d][XSTART];	\
      ny[d] = startny[d];			\
    }						\
  }						\
  y[tsxdim] = i_col = 0  			


#define STANDARD_ENDE_X				\
  FREE(xx)
 				

#define STANDARD_ENDE				\
  STANDARD_ENDE_X;				\
  FREE(yy)


#define GRIDCYCLES(SOME_FCTN)				\
  while (true) {					\
    SOME_FCTN;						\
    if (ygiven) { STANDARDINKREMENT_Y; RECYCLE_Y;}	\
    STANDARDINKREMENT_X;				\
  }					

#define GRIDCYCLE_X(SOME_FCTN)			\
  while (true) {				\
    SOME_FCTN;					\
    STANDARDINKREMENT_X;			\
  }
    
#define DO_INCREMENTY , yyy+=tsxdim, i_col++

#define PREPAREX							\
  MEMCOPY(x, xxx, xCpy);						\
  x[tsxdim] = i_row

#define PREPAREY				\
  PREPAREX;						\
  if (yyy < yend) { } else { yyy = y0; i_col = 0; }	\
  MEMCOPY(y, yyy, xCpy);					\
  y[tsxdim] = i_col

#define EMPTY 

#define NONGRIDCYCLE(INCREMENT, PREPARE, FCTN1, FCTN2)			\
  if (vdimSq == 1) {							\
    for (; i_row<tot; i_row++, xxx+=tsxdim INCREMENT){			\
      PREPARE;								\
      FCTN1;								\
    }									\
  } else {								\
    for (; i_row<tot; i_row++, xxx+=tsxdim INCREMENT){			\
      PREPARE;								\
      FCTN2;								\
    }									\
  }

#define PERFORM(UNIVAR_FCTN_X, MULTIVAR_FCTN_X, UNIVAR_FCTN, MULTIVAR_FCTN)  \
  if (grid) {					\
    if (ygiven || kernel) {						\
      STANDARDSTART_Y_SUPPL;						\
      if (vdimSq == 1) { GRIDCYCLES(UNIVAR_FCTN; value+=vdimSq);	\
      } else { GRIDCYCLES(MULTIVAR_FCTN); }				\
    } else { /* grid, y not given */					\
      if (vdimSq == 1) {						\
	GRIDCYCLE_X(UNIVAR_FCTN_X; value+=vdimSq);			\
      } else {GRIDCYCLE_X(MULTIVAR_FCTN_X); }				\
    }									\
  } else { /* not a grid */						\
    int localdim = tsxdim;						\
    double *xxx, *yyy = zero;						\
    if (trafo) {							\
      localdim = TransformLoc(cov, &xx, &yy, false);			\
      xxx = xx;								\
      if (ygiven) yyy = yy;						\
    } else {								\
      xxx=loc->x;							\
      if (ygiven) yyy=loc->y;						\
    }									\
    int /* xMem = sizeof(double) * (tsxdim + 1), */			\
      xCpy = sizeof(double) * localdim;					\
    assert(ygiven xor (yyy==zero));					\
    if (ygiven || kernel) {						\
      double *y0 = yyy,							\
	*yend = ygiven ? yyy + tsxdim * loc->ly : yyy;			\
      NONGRIDCYCLE(DO_INCREMENTY, PREPAREY, UNIVAR_FCTN; value+=vdimSq,	\
		   MULTIVAR_FCTN);					\
    } else {								\
      NONGRIDCYCLE(EMPTY, PREPAREX, UNIVAR_FCTN_X; value+=vdimSq,	\
		   MULTIVAR_FCTN_X);					\
    }									\
    if (err != NOERROR) XERR(err);					\
  }



void CovVario(model *cov, bool is_cov, bool pseudo, 
	      double *value);
///void CovIntern(int reg, double *x, double *value);
void CovIntern(int reg, double *x, double *y, Long lx, Long ly, double *value);
void PseudovariogramIntern(int reg, double *x, double *y,
			   Long lx, Long ly, double *value);
void PseudovariogramIntern(int reg, double *x, double *value);
void CovarianceMatrix(model *cov,double *v);
void CovarianceMatrixCol(model *Cov, int column, double *v);
void InverseCovMatrix(model *cov, double *v, double *det);
void StandardCovariance(model *cov, double *v);
void StandardCovMatrix(model *cov, double *v);
void StandardInverseCovMatrix(model *cov, double *inverse, double *det);
void StandardVariogram(model *cov, double *v);
void StandardPseudoVariogram(model *cov, double *v);


#endif

