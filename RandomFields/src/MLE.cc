/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 library for simulation of random fields 

 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
RO
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Linpack.h>
#include <Rmath.h>  
#include <stdio.h>  
//#include <stdlib.h>
#include <string.h>

#include "questions.h"
#include "primitive.h"
#include "operator.h"
#include "Processes.h"
#include "shape.h"
#include "rf_interfaces.h"


bool is_top(model *cov) {
  if (cov == NULL) BUG;
  return equalsnowInterface(cov) || isnowProcess(cov);
}



int required(double v, double *values, int n) {
  if (ISNA(v)) {
    for (int i=0; i<n; i++) if (ISNA(values[i])) return i;
  } else if (ISNAN(v)) { // in der Reihenfolge !!
    for (int i=0; i<n; i++) if (R_IsNaN(values[i])) return i;
  } else {
    for (int i=0; i<n; i++) if (!ISNA(values[i]) && v == values[i]) return i;
  }
  return MISMATCH;
}



int GetNAPosition(model *cov, double *values, int nvalues,
		  double *fillvalues, int nfillvalues, 
		   int *NAs, naptr_type mem,
		   //		   int *elmnts, elptr_type mem_elmnts,
		   NAname_type names, sortsofparam *sorts, 
		   int *vdim1, int *vdim2, int *found, bool *bayesian,
		  int * coord,
		   covptr_type covModels, 
		   int *covzaehler, int allowforintegerNA,
		   int SHORTlen, int printing, int depth, bool no_variance,
		   bool excludetrends
		   ) {
  /*
    printing <= 0 nothing
    1 only NA, NaN
    2 others as "*"
    3 others as value

    mem[0..NAs] : addresse wo (neue) NA position im Speicher	   
    sorts[0..NAs]     : see sortsofparam in RF.h
    found[0..NAs] : NA or NaN or any other value requested
  */

  // if (SHORTlen < 2) BUG;
  // printf("%s\n", NAME(cov));

  //  printing = 1;
  
  depth++;
  if (SHORTlen >= 255) SHORTlen = 254;
  KEY_type *KT = cov->base;
  SPRINTF(KT->error_loc, "'%.50s'", NICK(cov));

  model *sub;
  int i, c, r,
    constkappanr = -1,
    nr = COVNR,
    iso = OWNISO(0),
    namenr = 1,     
    *nrow =  cov->nrow,
    *ncol = cov->ncol;
  char shortname[255], shortD[255];
  bool Const = isConst(cov) && !isPlus(cov->calling);
  defn *C = DefList + nr,
    *CC = C; // nicht gatternr
  while(STRNCMP(C->name, InternalName, STRLEN(InternalName)) == 0) {
    C--;
    nr--;
  }
	 

  // C is pointer to true model; CC points to the first model passed
  // in initNerror.cc by IncludeModel or IncludePrim
  if (Const) {
    nr = CALLINGNR;
    CC = DefList + nr;
    model *calling = cov->calling;
    int kappas = CC->kappas;
    constkappanr = 0;
    while (constkappanr < kappas && calling->kappasub[constkappanr] != cov)
      constkappanr++;
    if (constkappanr < kappas) {
      while(STRNCMP(CC->name, InternalName, STRLEN(InternalName)) == 0) {
	CC--;
	nr--;
      }
    } else {
      Const = false;
      nr = COVNR;
      CC = C;
      covzaehler[nr]++; 
    }
  } else covzaehler[nr]++;
   
  strcopyN(shortname, CC->nick + 2, SHORTlen);
  if (covzaehler[nr] >= 2) {
    char dummy[255];
    strcopyN(dummy, shortname, SHORTlen-1);
    SPRINTF(shortname, "%.50s%d", dummy, covzaehler[nr]);
  }
  if (printing>0) { PRINTF("%s\n", CC->name); }
  // CC needed below for the kappa.names which are given

  SEXPTYPE *type = C->kappatype;
  for (i=0; i<C->kappas; i++) {
    bool istrend = equalsnowTrendParam(cov, i);
    // printf("%s %d %d %d\n", NAME(cov), i, istrend, excludetrends);
    if (istrend && excludetrends) continue;
    /*
    if (i==0 && type[i] == INTSXP && STRCMP(C->kappanames[i], ELEMENT) == 0){
      if (*elmnts >= MAX_MLE_ELMNTS ) 
	ERR("maximum number of models with variable elements reached");
      if (P0INT(i) == NA_INTEGER) ERR1("'%.50s' may not be NA", ELEMENT);
      mem_elmnts[(*elmnts)++] = PINT(i);
      // to do!! covModels[*NAs] hierher uebertragen !! s.u.

      continue;
    }
    */
    
    if ((sub = cov->kappasub[i]) == NULL || isnowRandom(sub)) {    
      if (nrow[i] == 0 || ncol[i] == 0) continue;
  
      if (sub != NULL) {
	int size = nrow[i] * ncol[i];
	for (r =0; r < size; P(i)[r++] = RF_NA);
      }

      if (printing > 0) {
	leer(depth); PRINTF("%s\n", C->kappanames[i]);
      }

      for (c=0; c<ncol[i]; c++) {
 	int nv = 0; // anzahl NA in aktuellem parameter
	for (r=0; r<nrow[i]; r++) {
	  if (*NAs >= MAX_NA) SERR("maximum number of NA reached");	
	  sorts[*NAs] = SortOf(cov, i, r, c, original_model);  // ANYPARAM;

	  if (isNotEstimable(sorts[*NAs])) continue;
	  
	  double v = 0.0; // dummy value
	  int idx = c * nrow[i] + r;

	  found[*NAs] = UNSET;
	  vdim1[*NAs] = r + 1;
	  vdim2[*NAs] = c + 1;
	  bayesian[*NAs] = sub != NULL;
	  coord[*NAs] = iso;
	
	  switch(type[i]) {
	  case REALSXP :
	    v = P(i)[idx];
	    mem[*NAs] = P(i) + idx;
	    found[*NAs] = required(v, values, nvalues);
	    // printf("%s %s %d %d %d NAs=%d %d v=%f\n", C->name, C->kappanames[i], r, c, found[*NAs], *NAs, nfillvalues, v);
	    if (*NAs < nfillvalues && found[*NAs]>=0)
	      mem[*NAs][0] = fillvalues[*NAs];
	    covModels[*NAs] = cov;
	    break;
	  case INTSXP : {
	    v = PINT(i)[idx] == NA_INTEGER ? RF_NA : (double) PINT(i)[idx];
	    if (allowforintegerNA) {
	      found[*NAs] = required(v, values, nvalues);
	    } else {
	      sortsofparam s = SortOf(cov, i, r, c, mle_conform);
	      if (!isNotEstimable(s) && required(v, values, nvalues) >= 0)
		SERR1("%20s: integer variables currently not allowed to be 'NA'", KNAME(i)); // !!!
	    }
	  }    
	    break;
	  case LISTOF + REALSXP : {
	    listoftype *q = PLIST(i);
	    double *p = q->lpx[r];
	    int j, 
	      end = q->nrow[r] * q->ncol[r];
	    for (j=0; j<end; j++)
	      if (ISNAN(p[j]))
		SERR("no NAs allowed in arguments that can be lists");
	    break;
	  }
	  default :
	    if (!isRObject(type[i]) && type[i]!=STRSXP) BUG;
	  }
	  
	  //	  isnan[*NAs] = R_IsNaN(v);
	  //	  if (ISNAN(v)) { // entgegen Arith.h gibt ISNA nur NA an !!

	  //	  printf("found = %d *NAs=%d\n", found[*NAs], *NAs);
	  
	  if (found[*NAs] >= 0) { // entgegen Arith.h gibt ISNA nur NA an !!
	    
	    if (printing > 0) {
	      if (nv>1 || (c>0 && printing > 1)) PRINTF(", ");
	      else leer(depth+1);
	    }
	    if (isDollar(cov)) {
	      // shortD partial name for R level
	      model *next = cov->sub[0];
	      while(isDollar(next) || isNatsc(next)) next = next->sub[0];// ok
	      if (covzaehler[NEXTNR] == 0) { // next wurde noch nicht
		// untersucht, somit ist covzaehler[NEXTNR] um 1 niedriger
		// als covzaehler[nr] !
		strcopyN(shortD, NICK(next) + 2, SHORTlen);
	      } else {
		char dummy[255];
		strcopyN(dummy, NICK(next) + 2, SHORTlen-1);
		SPRINTF(shortD, "%.50s%d", dummy, covzaehler[NEXTNR]+1);
	      }
	      if (i==DVAR) { 		
		vdim1[*NAs] = vdim2[*NAs] = r + 1;
		model *calling = cov->calling;
		if (calling == NULL) BUG;
		while (!is_top(calling) && (MODELNR(calling) ==PLUS))
		  calling = calling->calling;	
		if (calling == NULL) BUG;
		sorts[*NAs] = is_top(calling)
		  ? (equalsNugget(NEXTNR) ? NUGGETVAR : VARPARAM)
		  : ANYPARAM; 
		if (no_variance && sorts[*NAs] <= SIGNEDSDPARAM)
		  sorts[*NAs] = ANYPARAM;
		SPRINTF(names[*NAs], 
			sorts[*NAs] <= SIGNEDSDPARAM || sorts[*NAs] ==NUGGETVAR 
			? "%.50s.var" : "%.50s.vx",
			shortD); // for R level only
	      } else if (i==DSCALE) {
		sorts[*NAs] = SCALEPARAM;	
		SPRINTF(names[*NAs], "%.50s.s", shortD);// for R level only
	      } else {
		assert(i == DANISO);
		if (cov->q != NULL && cov->q[0] != 0.0) {
		  SERR("naturalscaling only allowed for isotropic setting");
		}
		sorts[*NAs] = r==c ? DIAGPARAM : ANISOPARAM;
		SPRINTF(names[*NAs], "%.50s.A.%d.%d", shortD, r + 1, c + 1);// R
	      }
	    } else { // not $
	      //printf("entering non $\n");
	      // standard setting !
	      // printf(">> %s %d %d\n", NAME(cov), i, equalsnowTrendParam(cov, i));
	      if (istrend) {// || nr==MLEMI XEDEFFECT)
		// printf("Trendparam\n");
		assert(type[i] == REALSXP);	      
		if (r == 0 && c==0) {
		  int rr, 
		    rowcol = nrow[i] * ncol[i];
		  for (rr=0; r<rowcol; r++) {
		    double vv = P(i)[rr];
		    if (required(vv, values, nvalues) < 0) 
		      SERR("in case of trend parameters either none or all values must be 'NA'");
		  }
		}

		//	printf("exlucde trend = %d\n", excludetrends);		
		sorts[*NAs] = TRENDPARAM;
	      } else {  	
		if (type[i] == INTSXP)
		  sorts[*NAs] = INTEGERPARAM; 
		else {
		  // r)ow, c)olumn; i:kappa, nr, C		  
		  //		print("%.50s k=%d no_var=%d, %d \n", C->name, i,
		  //		       no_variance, sorts[*NAs]);  
		  if (no_variance && sorts[*NAs] <= SIGNEDSDPARAM)
		    sorts[*NAs] = ANYPARAM;
		}
	      }
	      char kappashort[255];
	      if (isConst(cov) && !PisNULL(CONST_NAME)) {
		//		printf("******* *%s\n", P0CHAR(CONST_NAME));
		strcopyN(names[*NAs], P0CHAR(CONST_NAME),  2 * SHORTlen);
		namenr=1;
	      } else {
		strcopyN(kappashort, OWNKAPPA(CC, Const ? constkappanr : i),
			 SHORTlen);
		if (nr == DECLARE) SPRINTF(names[*NAs], ".%.50s", kappashort);
		else SPRINTF(names[*NAs], "%.50s.%.50s", shortname, kappashort);

		if (*NAs > 0 && 0==
		    STRNCMP(names[*NAs], names[*NAs-1], STRLEN(names[*NAs]))){
		  if (namenr == 1) {
		    SPRINTF(names[*NAs-1], "%.50s.%.50s.%d",
			    shortname, kappashort, namenr);
		  } 
		  namenr++; 
		  SPRINTF(names[*NAs], "%.50s.%.50s.%d", shortname,
			  kappashort,namenr); 
		} else namenr=1;
	      }
	    }
	  
	    if (printing > 0) {
	      if (printing <= 2) PRINTF("%d",*NAs + 1); 
	      else PRINTF("!%d", *NAs + 1);
	    }
	    //	    printf("%s %d\n", NAME(cov), *NAs);
	    (*NAs)++;
	    nv++;
	    // printf("NA's = %d\n", *NAs);
	  } else { // not is.na(v)
	    // kein NA an der Position; somit nur Anzeige
	    if (printing > 1) {
	      if (c>0) PRINTF(", "); else leer(depth+1);
	      if (printing <= 2)  PRINTF("*");
	      else {
		if (type[i]== REALSXP) PRINTF("%6.2e", v);
		else  PRINTF("%5d", (int) v);
	      }
	    }
	  }  // else kein NA 
	} //c
	if ((printing > 0 && nv > 0) || (printing > 1 && ncol[i] > 0)) {
	  PRINTF("\n");
	}
      } // r  
    } // sub == NULL || isnowRandom(sub); KEIN ELSE !!

    if (sub != NULL) {
      if (printing > 0) {
	leer(depth); PRINTF("%s\n", C->kappanames[i]);
      }
 
      sortsofparam sort = SortOf(cov, -i-1, 0, 0, original_model);
      int err = GetNAPosition(sub, values, nvalues, fillvalues, nfillvalues, 
			      NAs, mem,// elmnts, mem_elmnts,
			      names, sorts, vdim1, vdim2, found, bayesian,
			      coord, covModels,
			      covzaehler, allowforintegerNA, 
			      SHORTlen, printing, depth, 
			      sort != VARPARAM && sort != VARONLYMLE,
			      excludetrends);      
      if (err != NOERROR) RETURN_ERR(err);
      
      continue;
    }
  } // i
 
  for (i=0; i<MAXSUB; i++) {
    sub = cov->sub[i];

  
    if (sub != NULL) {
      if (printing > 0) {
	leer(depth);
	PRINTF("%s = ", C->subnames[i]);
      }
      int err = GetNAPosition(sub, values, nvalues, fillvalues, nfillvalues, 
			      NAs, mem, //elmnts, mem_elmnts,
			      names, sorts, vdim1, vdim2, found, bayesian,
			      coord, covModels,
			      covzaehler, allowforintegerNA, 
			      SHORTlen, printing, depth, 
			      false, excludetrends);
      if (err != NOERROR) RETURN_ERR(err);
  
    }
  }

  RETURN_NOERROR;

}


 
SEXP GetNAPositions(SEXP Model_reg, SEXP Model, SEXP x,
		    SEXP values, SEXP fillvalues,
		    // SEXP spatialdim, SEXP Time, SEXP xdimOZ,
		    SEXP integerNA, SEXP vdim, SEXP Print) {
 int Reg = INTEGER(Model_reg)[0];
 set_currentRegister(Reg);
 int i, //elmnts, 
    vdim1[MAX_NA], vdim2[MAX_NA], covzaehler[MAXNRCOVFCTS], found[MAX_NA],
     coord[MAX_NA];
  naptr_type mem;
  //  elptr_type mem_elmnts;
  covptr_type covModels;
  sortsofparam sorts[MAX_NA];
  bool bayesian[MAX_NA];  
  NAname_type names;
  SEXP ans;
  KEY_type *KT = KEYT();
  model *cov = KT->KEY[Reg];

  if (Model != R_NilValue) {
    bool skipchecks = GLOBAL_UTILS->basic.skipchecks;
    GLOBAL_UTILS->basic.skipchecks = true;
    cov = InitIntern(Reg, Model, x, true);
    GLOBAL_UTILS->basic.skipchecks = skipchecks;
  }
  KT = cov->base;

  SPRINTF(KT->error_loc, "getting positions with NA");

  
  if (length(values) == 0) {
    PROTECT (ans =  allocVector(INTSXP, 0)); 
  } else {
    int NAs = 0;
    for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);
   int err = GetNAPosition(cov,
			    REAL(values), length(values),
			    REAL(fillvalues), length(fillvalues),
			    &NAs, mem,//&elmnts, mem_elmnts,
			    names, sorts, vdim1, vdim2, found, bayesian, 
			    coord, covModels,
			    covzaehler,
			    INTEGER(integerNA)[0],
			    GLOBAL.fit.lengthshortname, 
			    INTEGER(Print)[0], 
			    0, false, true);
   SPRINTF(KT->error_loc, "'%.50s'", NICK(cov));
   OnErrorStop(err, cov);
 
    //    printf("found %d NAs\n", NAs);

    PROTECT (ans =  allocVector(INTSXP, NAs));
    for (i=0; i<NAs; i++) 
      INTEGER(ans)[i] = found[i]>=0 ? found[i] + 1 : NA_INTEGER;
  }
  UNPROTECT(1);
  // PMI0(cov);
  INTEGER(vdim)[0] = cov->vdim[0];
  //  printf("length values = %d %d\n", length(values), NA_INTEGER);
  return ans;
}



int countnas(model *cov, bool excludetrend, int level, sort_origin origin) {
    int i, r,  count, 
      *nrow =  cov->nrow,
      *ncol = cov->ncol;
    defn *C = DefList + COVNR; // nicht gatternr
  SEXPTYPE *type = C->kappatype;

  count= 0;
  for (i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] != NULL) {
      count += countnas(cov->kappasub[i],
			false,//trend can't appear anymore,so nothing to exclude
			level + 1, origin);
    }

    if (excludetrend &&	equalsnowTrendParam(cov, i)
	//&& !PisNULL(i) // wichtig, da ja noch der Parameter durch ein Modell 
	)
      // gegeben sein koennte
      continue;
    int pt = SortOf(cov, i, 0, 0, origin),
      endfor = nrow[i] * ncol[i];
    
    if (endfor == 0 || isNotEstimable(pt)) continue;
  
    if (type[i] == REALSXP) { 
     double *p = P(i);
      for (r=0; r<endfor; r++) if (ISNAN(p[r])) count++;
    } else if (type[i] == INTSXP) {
      int *p = PINT(i);
      for (r=0; r<endfor; r++) if (p[r] == NA_INTEGER) count++;
    } else {
      continue; // no NAs allowed
    }
  } // for
  
  model *sub;
  excludetrend &= isPlus(cov) ||
    (CALLING != NULL && CALLINGNR == PLUS && isMal(cov));
  for (i=0; i<MAXSUB; i++) {
    sub = cov->sub[i];
    if (sub == NULL) continue; 
    count += countnas(sub, excludetrend, level + 1, origin);
  }

  //  print("<< %d %d\n", level, count);

  return count;
}



void Take21internal(model *cov, model *cov_bound,
		    double **bounds_pointer, int *NBOUNDS) {
  /* determine the bounds for the parameters that 
     are to be estimated
   */

  int i, c, r,
    nv=0, // zaehlt NAs
    *nrow = cov->nrow,
    *ncol = cov->ncol,
    *nrow2 = cov_bound->nrow,
    *ncol2 = cov_bound->ncol;
  defn *C = DefList + COVNR; // nicht gatternr
  SEXPTYPE *type = C->kappatype;
  
  if (STRCMP(Nick(cov), Nick(cov_bound)) != 0) {
    ERR("models do not match.");
  }  

  for (i=0; i<C->kappas; i++) {
    if (cov->kappasub[i] != NULL) {
      Take21internal(cov->kappasub[i], cov_bound->kappasub[i],
		     bounds_pointer, NBOUNDS);
      continue;
    }

    if (equalsnowTrendParam(cov, i)) continue;
   
    int pt = SortOf(cov, i, 0, 0, original_model);
    if (C->kappatype[i] >= LISTOF || isNotEstimable(pt))  continue;
    
    if (nrow[i] != nrow2[i] || ncol[i] != ncol2[i]) {
        PRINTF("%s i: %d, nrow1=%d, nrow2=%d, ncol1=%d, ncol2=%d\n", 
	       C->name, i, nrow[i], nrow2[i], ncol[i], ncol2[i]);
  	ERR("lower/upper/user does not fit the model (size of matrix)");
    }
    
    for (r=0; r<nrow[i]; r++) {
      for (c=0; c<ncol[i]; c++) {
	  double v=RF_NA,
	      w=RF_NA; // value in aktuellem parameter
	int idx = c * nrow[i] + r;
	
// print("%d %d idx=%d %d %d %d nv=%d %d\n", r,c, idx, type[i], REALSXP, INTSXP,
//	      nv, *NBOUNDS);

	if (type[i] == REALSXP) {
	  v = P(i)[idx];
	  w = PARAM(cov_bound, i)[idx];
	}
	else if (type[i] == INTSXP) {
	  v = PINT(i)[idx] == NA_INTEGER ? RF_NA : (double) PINT(i)[idx];
	  w = PARAMINT(cov_bound, i)[idx] == NA_INTEGER 
	    ? RF_NA : (double) PARAMINT(cov_bound, i)[idx];	  
	}
	     

	if (ISNA(v)) { // ISNA gibt nur NA an !!
	  if (!isDollar(cov) || 
	      i == DVAR ||
	      (i== DSCALE && cov->q == NULL) || // ! natscaling 
	      i == DANISO		
	      ) // aniso ?? ABfrage OK ??
	   {
	     if (nv >= *NBOUNDS) {
	       PRINTF("%s %s, r=%d, c=%d: %d >= %d\n",
	       		C->name, C->kappanames[i], r, c, nv, *NBOUNDS);
  	       ERR("lower/upper/user does not fit the model (number parameters)");
	     }
	     (*bounds_pointer)[nv] = w;
	     nv++;
           }
	}  //ISNA
      } //c
    } // r
  } // i
  *NBOUNDS = *NBOUNDS - nv;
  (*bounds_pointer) += nv;

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      Take21internal(cov->sub[i], cov_bound->sub[i], bounds_pointer, NBOUNDS);
    }
  }
}


SEXP Take2ndAtNaOf1st(SEXP model_reg, SEXP Model, SEXP model_bound,
		      SEXP Totalpoints,
		      SEXP spatialdim, SEXP Time, SEXP XdimOZ, 
		      SEXP nbounds, SEXP skipchecks) {
  int m,
    d = 0,
    xdimOZ = INTEGER(XdimOZ)[0],
    xdimOZ3 =xdimOZ * 3,
    totalpoints = INTEGER(Totalpoints)[0],
    NBOUNDS_INT =INTEGER(nbounds)[0],
    *NBOUNDS= &NBOUNDS_INT,
    nr[2] = {INTEGER(model_reg)[0], MODEL_BOUNDS};
  if (nr[0] == nr[1]) RFERROR("do not use register 'model bounds'");
   
  SEXP bounds,
    models[2] = {Model, model_bound};
  bool oldskipchecks = GLOBAL_UTILS->basic.skipchecks;
  KEY_type *KT = KEYT();

#define MAXDIM0 4  
  double *bounds_pointer,
    T[3] = {0.0, 1.0, 1.0},
    X0[MAXDIM0 * 3],
    *X = xdimOZ <= MAXDIM0 ? X0 : (double*) MALLOC(xdimOZ3 * sizeof(double));
  // falls CheckModel abbricht, sammelt sich *X als Muell an
  // MAXDIM0 = 4 hat alle Standardfaelle abgedeckt
 
  X[d++] = 0.0;
  X[d++] = 1.0;
  X[d++] =  (double) totalpoints;
  while(d<xdimOZ3) {
    X[d++] = 0.0;
    X[d++] = 1.0;
    X[d++] = 1.0;
  }


 setKT_NAOK_RANGE(KT, true);
  if (LOGICAL(skipchecks)[0]) GLOBAL_UTILS->basic.skipchecks = true;
  for (m=1; m>=0; m--) { // 2er Schleife !!
    // assert(m==1);
    double *zero = ZERO(xdimOZ + LOGICAL(Time)[0], KT);

    //    printf("T3 = %f\n", T[2]);
    
    CheckModel(models[m], X, zero, T,
	       INTEGER(spatialdim)[0], // spatialdim
	       xdimOZ,
	       3, 0, // lx, ly
	       true, // grid
	       false, // distances
	       LOGICAL(Time)[0], 
	       R_NilValue, // xlist
	       KT, nr[m]);  
    //    TREE(key[nr[m]]);

    //   printf("XT3 = %f\n", T[2]);
 
    
    GLOBAL_UTILS->basic.skipchecks = oldskipchecks;
  }


  PROTECT(bounds = allocVector(REALSXP, *NBOUNDS));
  bounds_pointer = NUMERIC_POINTER(bounds);

  //  PMI(key[nr[0]]);
  //  PMI(key[nr[1]]);

  model **key = KT->KEY;
  Take21internal(key[nr[0]], key[nr[1]], &bounds_pointer, NBOUNDS);
  if (xdimOZ > MAXDIM0) FREE(X);

  if (*NBOUNDS != 0) RFERROR("lower/upper does not fit to model");
  UNPROTECT(1);
  return bounds;
}




void GetNARanges(model *cov, model *min, model *max, 
		 double *minpile, double *maxpile, int *NAs,
		 bool dosimulations, sort_origin origin) {
    /* 
       determine the ranges of the parameters to be estimated 
     */

  int i,  r,
    *nrow = cov->nrow,
    *ncol = cov->ncol;
  defn 
    *C = DefList + COVNR; // nicht gatternr
  SEXPTYPE *type = C->kappatype;
  double v, dmin, dmax;
  
  for (i=0; i<C->kappas; i++) {
    model *sub = cov->kappasub[i];

    int end = nrow[i] * ncol[i];    
  
    if (end > 0 && (sub == NULL || isnowRandom(sub))) {
      switch(type[i]) {
      case REALSXP :
	dmin = PARAM0(min, i);
	dmax = PARAM0(max, i);
	break;
      case INTSXP :
	dmin = PARAM0INT(min, i) == NA_INTEGER 
	  ? RF_NA : (double) PARAM0INT(min, i);
	dmax = PARAM0INT(max, i) == NA_INTEGER 
	  ? RF_NA : (double) PARAM0INT(max, i);
	break;
      case LISTOF + REALSXP : {
	int store = GLOBAL.general.set;
	GLOBAL.general.set = 0;
	dmin = LPARAM0(min, i);
	dmax = LPARAM0(max, i);
	GLOBAL.general.set = store;
	break;
      }
      default:
	if (isRObject(type[i]) || type[i] == STRSXP) {
	  dmin = 0.0;
	  dmax = 0.0;
	} else {
	  BUG;
	  dmin = dmax = RF_NA;
	}
      }

      if (sub != NULL && end == 1 && dosimulations) {// i.e. isnowRandom
	int simulations = 1000; 
	double rr,
	  minr = RF_INF,
	  maxr = RF_NEGINF;
	for (int k=0; k<simulations; k++) {
	  DORANDOM(sub, &rr);
	  if (minr > rr) minr = rr;
	  if (maxr < rr) maxr = rr;
	}
	if (minr > dmin) dmin = minr;
	if (maxr < dmax) dmax = maxr;
      }
      
      int pt = SortOf(cov, i, 0, 0, origin);
      if (isNotEstimable(pt) || equalsnowTrendParam(cov, i)) continue;
      
      for (r=0; r<end; r++) {
	v = RF_NA;
	if (type[i] == REALSXP) {
	  v = P(i)[r];
	} else if (type[i] == INTSXP) {
	  v = PINT(i)[r] == NA_INTEGER ? RF_NA : (double) PINT(i)[r];
	} else if (isRObject(type[i]) || type[i] == STRSXP ||
		   type[i] == LISTOF + REALSXP)
	  break; //  continue;  // !!!!!!!!!!!
	else BUG;
	
	if (ISNAN(v)) {
	  if (isDollar(cov)) {
	    assert(i!=DAUSER && i!=DPROJ);
	  }
	  minpile[*NAs] = dmin;
	  maxpile[*NAs] = dmax;
	  (*NAs)++;
	} // isna
      } // r
    } // sub == NULL | isnowRandom
    
    if (sub != NULL) {
      //  PMI(cov->kappasub[i]);
      GetNARanges(cov->kappasub[i], min->kappasub[i], max->kappasub[i], 
		  minpile, maxpile, NAs, dosimulations, origin);
    }
  } // kappas
 

  for (i=0; i<MAXSUB; i++) {
    if (cov->sub[i] != NULL) {
      GetNARanges(cov->sub[i], min->sub[i], max->sub[i], 
		  minpile, maxpile, NAs, dosimulations, origin);
    }
  }

  //  print("end getnaranges %d\n", *NAs);

}




sortsofeffect getTrendEffect(model *cov) {
  int nkappa = DefList[COVNR].kappas;
  for (int j=0; j<nkappa; j++) {
    if (equalsnowTrendParam(cov, j)) {
      if (!PisNULL(j)) {
	return ISNAN(P0(j)) ? FixedTrendEffect : DetTrendEffect;
      } else {
	if (cov->kappasub[j] == NULL) {	  
	  return DetTrendEffect; // 30.4.15 vorher FixedTrendEffect
	}
	else if (isnowRandom(cov->kappasub[j])) {
	  RFERROR("priors not allowed in the context of trends"); // return RandomEffect;
	} else if (COVNR == TREND && j == TREND_MEAN) return DetTrendEffect;
	else RFERROR("model too complex");
      }
    }
  }
  //  printf("here A \n");
  return DetTrendEffect;
}


int CheckEffect(model *cov) {
  int nr;    
  // bool na_var = false;
  double *p;
 
  if (COVNR == TREND) {
    int trend = TREND_MEAN, 
      effect = effect_error;
    bool isna ;
    
    if ((nr = cov->nrow[trend] * cov->ncol[trend]) > 0) {
      p = P(trend);
      isna = ISNAN(p[0]);
      if ((effect != effect_error) && ((effect==FixedTrendEffect) xor isna))
	SERR1("do not mix known effect with fixed effects in '%.50s'", 
	      NICK(cov));
      for (int i = 1; i<nr; i++) {
	if (ISNAN(p[i]) xor isna) 
	  SERR("mu and linear trend:  all coefficient must be known or all must be estimated");
      }
      effect = isna ? FixedTrendEffect : DetTrendEffect;
    } else if (cov->kappasub[trend] != NULL) effect = DetTrendEffect;
    return effect;
  }
    
      
  if (equalsnowTrend(cov)) {
    if (COVNR == MULT) {
      sortsofeffect current = getTrendEffect(cov->sub[0]);

      // printf("current %d\n", current);
      for (int i=1; i<cov->nsub; i++) {
	sortsofeffect dummy = getTrendEffect(cov->sub[i]);
	if (current != DetTrendEffect && dummy != DetTrendEffect)
	  ERR("trend parameter to be estimated given twice");
	if (current == DetTrendEffect) current = dummy;
      }
      if (current == effect_error) ERR("trend mismatch");
      return current;
    } else return getTrendEffect(cov);
  }
 
  return ErrorEffect;
  }


int GetEffect(model *cov,  likelihood_info *info, sort_origin origin) {
  if (isnowProcess(cov)) {
    assert(cov->key == NULL);
    assert(!PisNULL(GAUSS_BOXCOX));
    int nas = 0,
      total = cov->nrow[GAUSS_BOXCOX] * cov->ncol[GAUSS_BOXCOX];
    for (int i=0; i<total; i++)
      nas += ((bool) (ISNAN(P(GAUSS_BOXCOX)[i])));// see remark in Arith.h
    if (nas > 0) {
      if (info->neffect >= MAX_LIN_COMP) ERR("too many linear components");
      info->nas[info->neffect] = nas;
      info->effect[info->neffect] = DataEffect;
      (info->neffect)++;
    }
    return GetEffect(cov->sub[0], info, origin);
  }
  
  bool excludetrend = true,
    plus = COVNR == PLUS; // || COVNR == SELECT;
  int i,
    n = plus ? cov->nsub : 1;
  for (i=0; i<n; i++) {
    model *component = plus ? cov->sub[i] : cov;

    if (MODELNR(component) == PLUS) {
      GetEffect(component, info, origin);
      continue;
    }
    
    if (info->neffect >= MAX_LIN_COMP) ERR("too many linear components");
    info->effect[info->neffect] = CheckEffect(component);
    info->nas[info->neffect] = countnas(component, excludetrend, 0, origin);
    if (info->effect[info->neffect] == effect_error) 
      SERR("scaling parameter appears with constant matrix in the mixed effect part");
    if (info->effect[info->neffect] >= ErrorEffect || 
	info->effect[info->neffect] <= // 5.8.17 >=
	DataEffect) {
      //    printf("varmodel = %d %d %d  %d\n", info->varmodel,model_undefined, 
      //   model_morethan1, info->effect);
      info->varmodel =
	info->varmodel == model_undefined ? info->neffect : model_morethan1;
      info->Var = component;
    } 
    (info->neffect)++;
  }
 
  RETURN_NOERROR;
}
  



int SetAndGetModelInfo(model *key, int shortlen, 
		       int allowforintegerNA, bool excludetrend, // IN
		       int newxdim,  usr_bool globalvariance, // IN 
		       likelihood_info *info, sort_origin origin){
			// OUT
		       //  bool *trans_inv, bool *isotropic, double **Matrix,
		       //int *neffect, int effect[MAXSUB],
		       //int *NAs, int nas[MAXSUB], NAname_type names) {
  int i, rows,
    mem_vdim1[MAX_NA],
    mem_vdim2[MAX_NA],
    covzaehler[MAXNRCOVFCTS], 
    jump = UNSET,
    err = NOERROR;
  model 
    *cov =
    !equalsnowInterface(key) ? key : key->key == NULL ? key->sub[0] : key->key,
    *sub = !isnowProcess(cov) ? cov : cov->key == NULL ? cov->sub[0] : cov->key,
    *min=NULL,  *max=NULL, 
    *pmin=NULL,  *pmax=NULL, 
    *openmin=NULL,  *openmax=NULL;
  sortsofparam mem_sorts[MAX_NA];
  int mem_found[MAX_NA], mem_coord[MAX_NA];
  bool mem_bayesian[MAX_NA];
  double mle_min[MAX_NA], mle_max[MAX_NA], mle_pmin[MAX_NA], mle_pmax[MAX_NA],
    mle_openmin[MAX_NA], mle_openmax[MAX_NA],
    *matrix = NULL;
#define nvalues 2
  double values[nvalues] = { RF_NA, RF_NAN };
  mle_storage *s = key->calling != NULL ? key->calling->Smle : key->Smle;
  KEY_type *KT = cov->base; 
  assert(s != NULL);

  SPRINTF(KT->error_loc, "checking model");
 
  if (sub->pref[Nothing] == PREF_NONE) {
    err = ERRORINVALIDMODEL;
    goto ErrorHandling;
  }
  
  info->newxdim = newxdim;
  info->trans_inv = isXonly(SUB) && isnowNegDef(sub);
  info->isotropic = isAnyIsotropic(ISO(SUB, 0));
  if (info->trans_inv && info->isotropic) {
      info->newxdim = 1;
//	removeOnly(KEY() + modelnr); nicht wegnehmen
// da derzeit negative distanczen in GLOBAL.CovMatrixMult auftreten -- obsolete ?!
  }

  
  // GLOBAL_UTILS->basic.skipchecks = skipchecks;
  check_recursive_range(key, true);

  if ((err = get_ranges(cov, &min, &max, &pmin, &pmax, &openmin, &openmax))
      != NOERROR) goto ErrorHandling;

  //  s->ELMNTS = 
  s->NAS = 0;
  for (i=0; i<MAXNRCOVFCTS;  covzaehler[i++]=0);
  // !! ZWINGEND VOR  GetEffect DA BAYESSCHE VARIABLE AUF NA GESETZT WERDEN
  if ((err = GetNAPosition(cov, values, nvalues, NULL, 0,
			   &(s->NAS), s->MEMORY, 
		//s->ELMNTS , s->MEMORY_ELMNTS,
			   info->names, mem_sorts, mem_vdim1, mem_vdim2,
			   mem_found,
			   mem_bayesian,
			   mem_coord,
			   s->COVMODELS,
			   covzaehler,
			   allowforintegerNA,
			   shortlen,
			   (PL >= PL_COV_STRUCTURE) + (PL >= PL_DETAILS) + 
		           (PL >= PL_SUBDETAILS), 
			   0, false, excludetrend))
       != NOERROR) goto ErrorHandling;

  SPRINTF(KT->error_loc, "'%.50s'", NICK(key));

  //  PMI(sub);
  
  info->NAs = 0; GetNARanges(cov, min, max, mle_min, mle_max, &(info->NAs),
			     false, origin);
  //   PMI(cov);
  //   print("XXz modelnr=%d %d\n",info->NAs, s->NAS);
  if (info->NAs != s->NAS) BUG;
  info->NAs = 0; GetNARanges(cov, pmin, pmax, mle_pmin, mle_pmax, &(info->NAs),
			     true, origin);
  //print("XX-------- amodelnr=%d\n",info->NAs);
  if (info->NAs != s->NAS) BUG;
  info->NAs = 0; GetNARanges(cov, openmin, openmax, mle_openmin, mle_openmax,
			     &(info->NAs), false, origin);
  if (info->NAs != s->NAS) BUG;
 // print("XX-------- amodelnr=%d %d\n",info->NAs, s->NAS); BUG;

  assert(info->varmodel == model_undefined);
  assert(info->neffect == 0);
  // GetEffect ersetzt NA-variance durch 1.0 
  if ((err = GetEffect(cov, info, origin)) != NOERROR) goto ErrorHandling;
  if (info->NAs != s->NAS) BUG;
  //  

  //  print("XX-------- amodelnr=%d\n",info->NAs); BUG;
 
  rows = s->NAS; // bevor NA-variance durch 1.0 ersetzt wird
  info->globalvariance = false;
  info->pt_variance = NULL;
  assert(COVNR == GAUSSPROC || hasAnyEvaluationFrame(cov));

  if ((COVNR == GAUSSPROC || hasAnyEvaluationFrame(cov)) && 
      (globalvariance == Nan || globalvariance == True) &&
       info->varmodel != model_undefined &&
       info->varmodel != model_morethan1){
    // does a global variance exist?
     
    if (isDollar(info->Var) && !PARAMisNULL(info->Var, DVAR) &&
	info->Var->ncol[DVAR]==1 && info->Var->nrow[DVAR]==1) {
      double var = PARAM0(info->Var, DVAR);
      if ((info->globalvariance = (bool) (ISNAN(var)) &&
	   info->Var->kappasub[DVAR] == NULL)) {	
	info->pt_variance = PARAM(info->Var, DVAR);
	assert(info->pt_variance != NULL);
	*(info->pt_variance) = 1.0;
	(info->NAs)--;
	info->nas[info->varmodel]--;
	for (i=0; i<rows; i++) {
	  if (s->MEMORY[i] == info->pt_variance) {
	    jump = i;
	    for (int j=i+1; j<rows; j++) {
	      s->MEMORY[j-1] = s->MEMORY[j];
	      STRCPY(info->names[j-1],info->names[j]);
	    }
	    break;
	  }
	}
	if (i >= rows) BUG;
	s->NAS--;       
	rows--;
      }
    }
    info->globalvariance |= (globalvariance == true);    
  }
  if (info->NAs != s->NAS) BUG;
  s->PT_VARIANCE = info->pt_variance;

  if (rows == 0) {
    info->Matrix = NULL;
  } else {
    info->Matrix = (double*) MALLOC(rows * MINMAX_ENTRIES * sizeof(double));
    matrix = info->Matrix;
    for (int k=i=0; i<rows; i++, k++) {
      if (i == jump) k++;
      int j = i - rows;
      // printf("i=%d j=%d %10g %10g %d\n", i,j,  mle_pmin[i], mle_pmax[i],  mem_sorts[i]);1032
      
      matrix[j + MINMAX_PMIN * rows ] = mle_pmin[k];
      matrix[j + MINMAX_PMAX * rows] = mle_pmax[k];
      sortsofparam sort = mem_sorts[k];
      if (sort >= FIRSTONLYMLE && sort <= LASTONLYMLE) {
	switch (sort) {
	case VARONLYMLE : sort = VARPARAM; break;
	case CRITONLYMLE : sort = CRITICALPARAM; break;
	case ONLYMLE : sort = ANYPARAM; break;
	default : BUG;
	}
      }
      matrix[j + MINMAX_TYPE * rows] = sort;
      matrix[j + MINMAX_NAN * rows] = !ISNA(values[mem_found[k]]);
      matrix[j + MINMAX_MIN * rows] = mle_min[k];
      matrix[j + MINMAX_MAX * rows] = mle_max[k];
      matrix[j + MINMAX_OMIN * rows] = mle_openmin[k];
      matrix[j + MINMAX_OMAX * rows] = mle_openmax[k];
      matrix[j + MINMAX_ROWS * rows] = mem_vdim1[k];
      matrix[j + MINMAX_COLS * rows] = mem_vdim2[k];
      matrix[j + MINMAX_BAYES * rows] = mem_bayesian[k];
      matrix[j + MINMAX_COORD * rows] = mem_coord[k];
    }
  }

 ErrorHandling:
  
  if (min != NULL) COV_DELETE(&min, cov);
  if (max != NULL) COV_DELETE(&max, cov);
  if (pmin != NULL) COV_DELETE(&pmin, cov);
  if (pmax != NULL) COV_DELETE(&pmax, cov);
  if (openmin != NULL) COV_DELETE(&openmin, cov);
  if (openmax != NULL) COV_DELETE(&openmax, cov);
  
  RETURN_ERR(err);
}


SEXP SetAndGetModelInfo(SEXP Model_reg, SEXP Model, SEXP x, 
			int spatialdim, 
			bool distances, int lx, int ly, bool Time, int xdimOZ,
			int shortlen, int allowforintegerNA, 
			bool excludetrend,
			sort_origin origin) {
  //printf("set&get\n");
   // ex MLEGetNAPos
    /*
      basic set up of covariance model, determination of the NA's
      and the corresponding ranges.
     */
//  bool skipchecks = GLOBAL_UTILS->basic.skipchecks;
  int 
    // NAs,
    //   unp = 0,
    //  effect[MAXSUB], nas[MAXSUB],
    //*pnas = nas,
    // neffect = 0,
    err = NOERROR;
  // double *Matrix = NULL;
  const char *colnames[MINMAX_ENTRIES] =
    {"pmin", "pmax", "type", "NAN", "min", "max", "omin", "omax",
     "row", "col", "bayes", "iso"};
  SEXP 
   matrix, nameAns, nameMatrix, RownameMatrix, 
   ans=R_NilValue;
  bool 
    do_not_del_info, isList_x = TYPEOF(x) == VECSXP;
  likelihood_info local,
    *info = NULL;
 //  isListofList_x = isList_x;
 
  if (isList_x) {
    if (length(x) == 0) BUG;
   //  isListofList_x = TYPEOF(VECTOR_ELT(x, 0)) == VECSXP;
  }
 
  int Reg = INTEGER(Model_reg)[0];
  set_currentRegister(Reg);
  KEY_type *KT = KEYT();
  double *zero = ZERO(xdimOZ + Time, KT), *x0 = zero, *y0 = zero;
  setKT_NAOK_RANGE(KT, true);
  if (length(x) != 0) {
    x0 = isList_x ? NULL : REAL(x);
    y0 = NULL;
  }
  //  printf("set^get A\n");
  CheckModel(Model, x0, y0, y0,
	     spatialdim, xdimOZ,
	     lx, ly, // lx, ly
	     false, //grid
	     distances, //distances
	     Time, 
	     isList_x ? x : R_NilValue,
	     KT, Reg);
 
  //  printf("set^get B\n");
  model *cov = KT->KEY[Reg],
    *Likeli = cov;
  //TREE(cov);
 
  if (equalsnowInterface(Likeli)) {
    model *process = cov->key != NULL ? Likeli->key : Likeli->sub[0];
    if (Likeli->Slikelihood == NULL && isnowProcess(process))
      Likeli = process;
  }

  
  ONCE_NEW_STORAGE(mle);
  mle_storage *s = cov->Smle;

  
  //  PMI0(Likeli);
  //  printf("set^get C\n");

  if ((do_not_del_info = Likeli->Slikelihood != NULL)) {    
    info = &(Likeli->Slikelihood->info);
    // pmi(cov,1);   BUG;
  } else {
    info = &local;
    likelihood_info_NULL(info);
    err = SetAndGetModelInfo(cov, shortlen, allowforintegerNA,// no PROTECT( needed
			     excludetrend, 
			     xdimOZ, GLOBAL.fit.estimate_variance,
			     info, origin);
    OnErrorStop(err, cov);
  }
  assert(s != NULL);

  int rows, NaNs;
  double *infoNaN;
  rows = s->NAS;
  NaNs = 0;
  infoNaN = info->Matrix + rows * (MINMAX_NAN - 1); // MINMAX_NAN is R coded
  for (int i=0; i<rows; i++) {
    //printf("i=%d %10g %d; %d %d %d %d\n", i, infoNaN[i], rows,MINMAX_PMIN,  MINMAX_PMAX, MINMAX_TYPE, MINMAX_NAN);
    NaNs += infoNaN[i];
  }
  //  printf("set^get CD\n");
  assert(s != NULL);
  
  PROTECT(matrix =  allocMatrix(REALSXP, rows, MINMAX_ENTRIES));
  PROTECT(RownameMatrix = allocVector(STRSXP, rows));
  PROTECT(nameMatrix = allocVector(VECSXP, 2));
#define total 8
  PROTECT(ans = allocVector(VECSXP, total));
  PROTECT(nameAns = allocVector(STRSXP, total));

  //printf("%d\n", info->Matrix != NULL);

  if (rows > 0) {
    MEMCOPY(REAL(matrix), info->Matrix, rows * MINMAX_ENTRIES * sizeof(double));
    for (int i=0; i<rows; i++) {
      SET_STRING_ELT(RownameMatrix, i, mkChar(info->names[i]));
    }
  }
  //  printf("set^get CE\n");
  SET_VECTOR_ELT(nameMatrix, 0, RownameMatrix);
  SET_VECTOR_ELT(nameMatrix, 1, Char(colnames, MINMAX_ENTRIES));
  setAttrib(matrix, R_DimNamesSymbol, nameMatrix);
   
  int i;
  i = 0;
  SET_STRING_ELT(nameAns, i, mkChar("minmax"));
  SET_VECTOR_ELT(ans, i++, matrix);
  SET_STRING_ELT(nameAns, i, mkChar("trans.inv")); // !!
  SET_VECTOR_ELT(ans, i++, ScalarLogical(info->trans_inv));  
  SET_STRING_ELT(nameAns, i, mkChar("isotropic"));
  SET_VECTOR_ELT(ans, i++, ScalarLogical(info->isotropic));  
  SET_STRING_ELT(nameAns, i, mkChar("effect"));
  SET_VECTOR_ELT(ans, i++, Int(info->effect, info->neffect));  
  SET_STRING_ELT(nameAns, i, mkChar("NAs"));
  SET_VECTOR_ELT(ans, i++, Int(info->nas, info->neffect));  
  SET_STRING_ELT(nameAns, i, mkChar("NaNs"));
  SET_VECTOR_ELT(ans, i++, ScalarInteger(NaNs));  
  SET_STRING_ELT(nameAns, i, mkChar("xdimOZ"));
  SET_VECTOR_ELT(ans, i++, ScalarInteger(info->newxdim));  
  SET_STRING_ELT(nameAns, i, mkChar("matrix.indep.of.x"));
  SET_VECTOR_ELT(ans, i++,
		 ScalarLogical(cov->matrix_indep_of_x));  
  setAttrib(ans, R_NamesSymbol, nameAns);
  assert(i==total);

  UNPROTECT(5);
  if (info != NULL && !do_not_del_info) likelihood_info_DELETE(info);

  return ans;

}


 
 
SEXP SetAndGetModelLikelihood(SEXP Model_reg, SEXP Model, SEXP x, SEXP origin) {
  //  printf("set^get\n");
  return SetAndGetModelInfo(Model_reg, Model, x, 0,			    
			    false, 0, 0, false, 0,// only dummy variables
			    GLOBAL.fit.lengthshortname, LIKELI_NA_INTEGER,
			    LIKELI_EXCLUDE_TREND,
			    (sort_origin) INTEGER(origin)[0]);
}



SEXP SetAndGetModelInfo(SEXP Model_reg, SEXP Model, SEXP spatialdim, 
			SEXP distances, SEXP ygiven, SEXP Time,
			SEXP xdimOZ, SEXP shortlen, SEXP allowforintegerNA, 
			SEXP excludetrend) {
  return SetAndGetModelInfo(Model_reg, Model, R_NilValue, 
			    INTEGER(spatialdim)[0], LOGICAL(distances)[0], 1,
			    LOGICAL(ygiven)[0] ? 1 : 0,
			    LOGICAL(Time)[0], INTEGER(xdimOZ)[0], 
			    INTEGER(shortlen)[0],
			    INTEGER(allowforintegerNA)[0],
			    LOGICAL(excludetrend)[0],
			    original_model);
}





void expliciteDollarMLE(int* reg, double *values) { // 
    // userinterfaces.cc 
  model *key = KEY()[*reg];
  mle_storage *s = key->Smle; assert(s != NULL);
  int i, un,
    NAs = s->NAS;
 	
  // first get the naturalscaling values and devide the preceeding 
  // scale model by this value 
  if (GLOBAL.general.naturalscaling==NATSCALE_MLE) {
    iexplDollar(key, true);
  }

  // Then get the values out of the model
  for (un=i=0; un<NAs; i++) {
    values[un++] = s->MEMORY[i][0];
    s->MEMORY[i][0] = RF_NA;
  }
}



void PutValuesAtNAintern(int *reg, double *values, bool init){
  model *key = KEY()[*reg];
  mle_storage *s = key->Smle; assert(s != NULL);
  int i, un,
    NAs = s->NAS;
  defn *C = NULL;
  model *cov = NULL;
  gen_storage S;
  double *pt_variance = s->PT_VARIANCE;
  gen_NULL(&S);
  S.check = S.dosimulate = false;
  // set ordinary parameters for all (sub)models
  
  set_currentRegister(*reg);
  for (un=i=0; i<NAs; i++) {
 //      print("reg=%d i=%d %d %ld %10g %ld\n", *reg, i, NAs, s->MEMORY[i], values[un], pt_variance);
    //     print("mem=%ld %10g\n", (Long) s->MEMORY[i], values[un]) ;  
    if (s->MEMORY[i] == pt_variance) BUG; //continue;
    s->MEMORY[i][0] = values[un++];
  }
   
   //   PMI(KEY()[*reg]);
    
  if (init)
    for (i=0; i<NAs; i++) {
      cov = s->COVMODELS[i];
      C = DefList + COVNR;       
      if (i==0 || cov != s->COVMODELS[i-1]) {
	if (!isDummyInit(C->Init)) {
	  assert(cov->checked);
	  // if (!(cov->initialised || isnowVariogram(cov))) TREE0(cov);
	  assert(cov->initialised || isnowVariogram(cov));
	  C->Init(cov, &S);
	}
      }
      //  
    }

  //APMI(KEY()[*reg]);
  
  //printf("done\n");

  //  int one = 1;  setListElements(reg, &one, &one, &one);
}

void PutValuesAtNA(int *reg, double *values){
  PutValuesAtNAintern(reg, values, true);
}

void PutValuesAtNAnoInit(int *reg, double *values){
  PutValuesAtNAintern(reg, values, false);
}
