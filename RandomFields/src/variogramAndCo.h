

/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2017  Martin Schlather

This
 program is free software; you can redistribute it and/or
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

#define FINISH_START(whereSfctn, whereVdim, ignore_x, ignore_y, idxVdim1, \
		     TOTX, TOTY, I_COL_BASE)				\
  /*//  printf("ignor %d %d\n", ignore_x, ignore_y); */			\
  GETSTORAGE(fctn, whereSfctn, fctn);					\
  if (fctn == NULL) { /* // PMI0(whereSfctn);printf("\n\n\n\n");PMI0(whereVdim); printf("fctn=NULL\n"); */ BUG;} \
  assert(Loc(whereSfctn) != NULL);					\
  double *caniso = LocAniso(whereSfctn);				\
  bool grid = LocgridY(whereSfctn, !(ignore_x)) && caniso==NULL,	\
    Time = LocTime(whereSfctn),						\
    trafo = ((Time && !grid)|| caniso != NULL);				\
  grid &= !trafo;							\
  if (trafo && LocDist(whereSfctn))					\
    ERR("complex construction with distances currently not programed -- pls contract maintainer"); /* ODER EHER BUG WENNS AUFTRITT? TO DO */ \
   Long totX = (TOTX) ? (TOTX) : LoctotalpointsY(whereSfctn, !(ignore_x)), \
    totY = (TOTY) ? (TOTY) : LoctotalpointsY(whereSfctn, ignore_y);	\
  int d,								\
    i_col_base = I_COL_BASE,						\
    *start=fctn->start,							\
    *end=fctn->end,							\
    i_row = 0,								\
    *nx=fctn->nx,							\
    tsxdim = PREVTOTALXDIM,						\
    VARIABLE_IS_NOT_USED localdim = tsxdim; /* dim after trafo */	\
  if (end == NULL) {/* //printf("end=NULL\n");*/ BUG;}			\
  double *x = fctn->x,							\
    *xstart= fctn->xstart,						\
    *inc=fctn->inc,							\
    *xx = NULL;		/* dito			   */			\
  assert(x != NULL);							\
  coord_type gr = LocgrY(whereSfctn, !(ignore_x));			\
  if (grid) {								\
    STANDARDSTART_X;							\
  } else x = LocY(whereSfctn, !ignore_x);				\
  									\
  /*// printf("************************************************** %s %d %d %d\n", #whereSfctn, (ignore_x), (ignore_y), LocHasY(whereSfctn)); */ \
  bool yAvailable = ((ignore_x) || !(ignore_y)) && LocHasY(whereSfctn),/* might be changed !*/ \
    extraData =  !((ignore_x) xor (ignore_y)) && LocHasY(whereSfctn),	\
    xswapped = (ignore_x) && LocHasY(whereSfctn),			\
    ggridY = LocgridY(whereSfctn, ignore_y),		\
    trafoY = ((Time && !ggridY) || caniso != NULL);			\
  ggridY &= !trafoY;							\
  int i_col = 0,							\
    vdim0  = (whereVdim)->vdim[0],					\
    vdim1 = (whereVdim)->vdim[idxVdim1],				\
    vdimSq = vdim0 * vdim1, /*not nec. squ!!*/				\
    *endy =fctn->endy,							\
    *startny =fctn->startny,						\
    *ny = fctn->ny;		 					\
		 							\
  Long vdim0totX = (Long) vdim0 * totX,					\
    vdim1totX = vdim1 * totX,						\
    VARIABLE_IS_NOT_USED VDIM_0, VARIABLE_IS_NOT_USED NEND,		\
    VARIABLE_IS_NOT_USED  NINCR, VARIABLE_IS_NOT_USED  MINCR,		\
    VARIABLE_IS_NOT_USED ENDFORINCR,					\
    VARIABLE_IS_NOT_USED VDIMtotX,					\
    vdimSqtotX = vdimSq * totX,						\
    vdim0totSq = vdim0totX * totY,					\
    vdimSqtotSq = vdimSq * totX * totY,					\
    err = NOERROR;							\
  assert(fctn->y != NULL);						\
  double *y = fctn->y,							\
    *ystart =fctn->ystart,						\
    *yy = NULL,		/* set by TransformLoc; freed  */		\
    *incy =fctn->incy;							\
  coord_type grY = LocgrY(whereSfctn, ignore_y);			\
    if (global->general.vdim_close_together) {				\
    /* v-dimensions close together */					\
    VDIM_0 = vdim0;							\
    VDIMtotX = vdimSqtotX;						\
    NINCR = vdim0totX;							\
    MINCR = 1;								\
    NEND = vdimSqtotX;							\
    ENDFORINCR = vdim0;							\
  } else {								\
    /* values of any *single* multivariate component close together */	\
    /* default in global->CovMatrixMulti */				\
    VDIM_0 = 1;								\
    VDIMtotX = vdim0totX;						\
    NINCR = vdim0totSq;							\
    ENDFORINCR = vdim0totX;						\
    NEND = vdimSqtotSq;							\
    MINCR = totX;							\
  }									\
  /*  PMI1(whereSfctn);	 PMIR(whereVdim); */				\
  DEFAULT_INFO(info);							\
  info[INFO_IDX_X] = info[INFO_IDX_Y] = 0



#define STANDARDSTART_X				\
  assert(Loc(cov) != NULL && gr[0] != NULL);	\
  for (d=0; d<tsxdim; d++) {			\
    inc[d] = gr[d][XSTEP];			\
    start[d] = 0;				\
    end[d] = (int) gr[d][XLENGTH];		\
    nx[d] = start[d];				\
    x[d] = xstart[d] = gr[d][XSTART];		\
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
  info[INFO_IDX_X] = ++i_row;	\
  if (d < tsxdim) { } else break		
  // end define StandardInkrement
        

// ny wird nur inv CovMatrix verwendet!
#define STANDARDSTART_Y_SUPPL						\
  for (d=0; d<tsxdim; d++){						\
    incy[d] = grY[d][XSTEP];						\
    y[d] = ystart[d] = grY[d][XSTART];					\
    endy[d] = (int) grY[d][XLENGTH];					\
    ny[d] = startny[d] = 0;						\
  }									\
  i_col = 0;							\
  info[INFO_IDX_Y] = i_col + i_col_base


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
  i_col++;					\
  info[INFO_IDX_Y] = i_col  + i_col_base		     
   // end define Stand ardInkrement_Y


#define RECYCLE_Y					\
  if (d < tsxdim) { } else {				\
    for (d=0; d<tsxdim; d++) {				\
      y[d] = ystart[d];		\
      ny[d] = startny[d];				\
    }							\
    i_col = 0;					\
    info[INFO_IDX_Y] = i_col + i_col_base;	\
  }


#define STANDARD_ENDE_X				\
  FREE(xx)


#define STANDARD_ENDE				\
  STANDARD_ENDE_X;				\
  FREE(yy)



#define GRIDCYCLES(SOME_FCTN)				\
  while (true) {					\
    SOME_FCTN;						\
    if (totY > 1) { assert(yAvailable &&  ggridY); STANDARDINKREMENT_Y; RECYCLE_Y;}	\
    STANDARDINKREMENT_X;	\
  }				

#define GRIDCYCLE_X(SOME_FCTN)			\
   while (true) {				\
     SOME_FCTN;					\
    STANDARDINKREMENT_X;			\
  }



    
#define DO_INCREMENTY  if (totY > 1) { i_col++;  y+=tsxdim;}

#define PREPAREX							\
  info[INFO_IDX_X]  = i_row

#define PREPAREY				\
  PREPAREX;						\
  if (y < yend) { } else { y = y0; i_col = 0; }	\
  info[INFO_IDX_Y] = i_col + i_col_base;

#define EMPTY 

#define NONGRIDCYCLE(INCREMENT, PREPARE, FCTN1, FCTN2)			\
  if (vdimSq == 1) {							\
    for (; i_row<totX; i_row++, x+=tsxdim){				\
      PREPARE;								\
      FCTN1;								\
      INCREMENT								\
    }									\
  } else {								\
    for (; i_row<totX; i_row++, x+=tsxdim){				\
      PREPARE;								\
      FCTN2;								\
      INCREMENT								\
    }									\
  }

#define PERFORM_PREPARE							\
  if (trafo) { /* implies no grid */					\
    localdim = TransformLoc(cov, NULL, xswapped ? NULL : &xx,		\
			    xswapped ? &xx : NULL,  False);		\
    x = xx;								\
  }									\
  if (yAvailable && !ggridY){						\
    if (trafoY) { /* falls xswapped && yAvailable, so auch yswapped	*/	\
      localdim = TransformLoc(cov, NULL, xswapped ? &yy : NULL,\
			      xswapped ? NULL : &yy, False);		\
      y = yy;								\
    } else y=LocY(cov, xswapped);					\
  } else if (yAvailable && ggridY) {					\
    STANDARDSTART_Y_SUPPL;						\
  } else y=zero
	   
	   
#define PERFORM(UNIVAR_FCTN_X, MULTIVAR_FCTN_X, UNIVAR_KERNEL,		\
		MULTIVAR_KERNEL)					\
  assert(!yAvailable || totY==1 || !(grid xor ggridY));			\
  assert(!extraData || callingKernel);					\
  if (grid) {								\
     if (extraData || callingKernel) {						\
      if (vdimSq == 1) { GRIDCYCLES(UNIVAR_KERNEL; v++); }		\
      else { GRIDCYCLES(MULTIVAR_KERNEL); }				\
    } else { /* grid, y not given */					\
      if (vdimSq == 1) { GRIDCYCLE_X(UNIVAR_FCTN_X; v++); }		\
      else { GRIDCYCLE_X(MULTIVAR_FCTN_X); }				\
    }									\
  } else { /* x not a grid */						\
    assert(yAvailable xor (y==zero));					\
    if (extraData || callingKernel) {					\
      double *y0 = y,							\
	*yend = yAvailable ? y + localdim * totY : y;			\
      NONGRIDCYCLE(DO_INCREMENTY, PREPAREY,				\
		   UNIVAR_KERNEL; v+=vdimSq, MULTIVAR_KERNEL);		\
    } else {								\
      NONGRIDCYCLE(EMPTY, PREPAREX, UNIVAR_FCTN_X; v+=vdimSq,		\
		   MULTIVAR_FCTN_X);					\
    }									\
    if (err != NOERROR) XERR(err);					\
  }


// kein komma in der nachfolgenden Defintion!!
#define VDIM_LOOP(DO)							\
   double *C = v + VDIM_0 * i_row + i_col * VDIMtotX;			\
  int m = 0;								\
  for (Long n1=0; m<vdim1; n1+=NINCR) { /* spaltenschleife */		\
    /* n1: index fuer ergebnismatrix C fuer spaltenschleife */		\
    int u = m * vdim0;/*richtiger idx der cross matrix fuer zeilenschleife*/ \
    int w = m; /* index der (virtuell) transponierten cross Matrix */	\
    for (Long m1=n1; w<vdimSq; m1+=MINCR) { /* zeilenschleife */	\
      /* printf("%d %d; %d %d; row=%d, col%d\n", m1, VDIM_0 * i_row + i_col * VDIMtotX + m1, totX, totY,  i_row, i_col); assert( VDIM_0 * i_row + i_col * VDIMtotX + m1< 400); //*/ \
      C[m1] = DO;							\
      u++; w+=vdim1;							\
    }									\
    m++;								\
  }


void PseudovariogramIntern(int reg, double *x, double *y,
			   Long totalpoints, Long totalpointsY, double *v);
void PseudovariogramIntern(int reg, double *x, double *v);
void CovarianceMatrix(model *cov, bool ignore_y, double *v);
void CovarianceMatrix(model *Cov, bool ignore_y, int *idx, int Nidx,
		      double *v);
void CovarianceMatrixCols(model *Cov, bool ignore_y, int row, double *v);
void StandardCovMatrix(model *cov, bool ignore_y, double *v);
//void InverseCovMatrix(model *cov, double *v, double *det);
//void StandardInverseCovMatrix(model *cov, double *inverse, double *det);
void Variogram(model *cov, double *v);
void Covariance(model *cov, double *v);
void CovarianceT(model *cov, int base_i_col, double *v);
void Pseudomadogram(model *cov, double alpha, double *v);

#endif

