
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 (library for simulation of random fields)
 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redist ribute it and/or
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

#include <R.h>
#include <Rdefines.h>
#include <stdio.h>  
 #include <string.h>
#include <R_ext/Linpack.h>

#include "def.h"
#include "questions.h"
#include "primitive.h"
#include "Coordinate_systems.h"
#include "Processes.h"
#include "xport_import.h"
#include "operator.h"
#include "PoissonPolygon.h"
#include "families.h"


#define ALL_NULL(name)				\
  void name##_NULL(name##_storage *x){		\
    if (x == NULL) return;			\
    MEMSET(x, 0, sizeof(name##_storage));	\
  }						

#define ALL_NULL_BUT(name, exceptions)		\
  void name##_NULL(name##_storage *x){		\
    if (x == NULL) return;			\
    MEMSET(x, 0, sizeof(name##_storage));	\
    exceptions;					\
  }						



void LOC_SINGLE_NULL(location_type *loc, int len, int dim) {
  int d;
  loc->spatialdim = loc->timespacedim = loc->xdimOZ = loc->rawset = UNSET;
  loc->xgr = (double **) MALLOC(dim * sizeof(double *));  
  loc->grY = (double **) MALLOC(dim * sizeof(double *));
  for (d=0; d<dim; d++) {
    loc->xgr[d] = loc->grY[d] = NULL;
  }
  loc->totalpoints = loc->spatialtotalpoints = 0;
  loc->grid = loc->gridY = loc->distances = loc->Time = false;
  loc->delete_x = loc->delete_y = true;
  loc->x = loc->Y = loc->caniso = NULL;
  loc->rawidx = NULL;
  //  loc->T[0] = loc->T[1] = loc->T[2] = 0.0;
  //loc->i_row = loc->i_col = I_COL_NA;
  loc->cani_ncol = loc->cani_nrow = NA_INTEGER;
  LocLocSets(loc) = len;
}
     

void LOC_NULL(location_type **Loc, int len, int dim) {
  int i;
  for (i=0; i<len; i++) LOC_SINGLE_NULL(Loc[i], len, dim);
}
 


location_type **LOCLIST_CREATE(int n, int dim) {
  int i;
  location_type **Loc = NULL;

  Loc = (location_type**) CALLOC(n, sizeof(location_type*));
  
  //  printf("LOCLIST_CREATE(int n, %d %ld\n", n, Loc);

  for (i=0; i<n; i++) Loc[i] = (location_type*) MALLOC(sizeof(location_type));
  LOC_NULL(Loc, n, dim);
  assert(LocPSets(Loc) >= 1 && LocPSets(Loc) <= 1e6);
  return Loc;
}

void LOCLIST_CREATE(model *cov, int n) {
  cov->ownloc = LOCLIST_CREATE(n, OWNXDIM(0));
}

void LOC_SINGLE_DELETE_Y(location_type *loc) {
  // called hier and in getNset.cc
  assert(loc != NULL);
  if (loc != NULL) { // sicherheitshalber
    if (loc->totalpointsY > 0 && loc->gridY) {
      assert(loc->grY != NULL);
      FREE(loc->grY[0]);
      // do not  FREE(loc->grY); here!!
    }
    if (loc->delete_y) FREE(loc->Y);
    loc->totalpointsY = 0;
  }
}
  

void LOC_SINGLE_DELETE(location_type **Loc) {
  assert(Loc != NULL);
  //printf("LOCLIST_SINGLE_DELE %ld\n", Loc); 
  location_type *loc = *Loc;    
  if (loc != NULL) {
    LOC_SINGLE_DELETE_Y(loc); // very first to deleted!
     FREE(loc->caniso);
    FREE(loc->rawidx);
    // it may happen that both are set ! Especially after calling 
    // partial _loc _set in variogramAndCo.cc
    if (loc->grid) {
      assert(loc->xgr != NULL);
      FREE(loc->xgr[0]);
    }
    if (loc->delete_x) FREE(loc->x); 

    FREE(loc->xgr); // always to be freed!
    FREE(loc->grY); // always to be freed, but not in LOC_SINGLE_DELETE_Y

    UNCONDFREE(*Loc);
  }
}



void LOC_DELETE(location_type ***Loc) { // OK
  if (*Loc == NULL) return;
  //  printf("LOCLIST_DELE %ld\n", *Loc); 
  int i, 
    len = (*Loc)[0]->len;
  for (i=0; i<len; i++) {
    LOC_SINGLE_DELETE( (*Loc) + i );    
  }

  assert(**Loc == NULL); 
  UNCONDFREE(*Loc);
}


listoftype * LIST_CREATE(int len, int Rtype) {
  if (len <= 0) BUG;
  listoftype* q = (listoftype *) MALLOC(sizeof(listoftype));
  q->lpx = (double **) CALLOC(len, sizeof(double *));
  q->ncol = (int*) CALLOC(len, sizeof(int));
  q->nrow = (int*) CALLOC(len, sizeof(int));
  q->deletelist = true;
  q->len = len;
  q->Rtype = Rtype;
  return q;
}

  

void LIST_DELETE(listoftype **x) {
  if (x == NULL) return;
  listoftype *q = *x;
  if (q != NULL) {
    assert(q->lpx != NULL);
    if (q->deletelist) {
      for (int i=0; i<q->len; i++) {
	FREE(q->lpx[i]);
      }
      FREE(q->lpx);
      FREE(q->ncol);
      FREE(q->nrow);
    }
    FREE(*x);
  }
}



void MPPPROPERTIES_NULL(mpp_properties *Mpp) {
  // Mpp->refradius= 
  int i;
  for (i=0; i<MAXMPPVDIM; i++) Mpp->maxheights[i] = RF_INF; // maxv
  Mpp->unnormedmass = RF_NA;
  Mpp->mM = Mpp->mMplus = NULL;
  //   Mpp->refsd = 
  //Mpp->totalmass = RF_NA;
  //Mpp->methnr = UNSET;
  //Mpp->loc_done = false;
}

void MPPPROPERTIES_DELETE(mpp_properties *Mpp) {   
  FREE(Mpp->mM);
  FREE(Mpp->mMplus);
}

 

void COV_DELETE_WITHOUTSUB(model **Cov, model *save) {
  // DIESES DING IST GEFAEHLICH DA prevloc ff NICHT BEHANDELT IST
  // save : where to save the error message constained in Cov when Cov is deleted?
  model *cov = *Cov;
  assert(cov != NULL);
  //  printf("deleting %.50s\n", NAME(cov));
  
  int i, j,
    lastparam = (COVNR < 0) ? MAXPARAM : DefList[COVNR].kappas; 
  for (i=0; i<lastparam; i++) {
    int type = DefList[COVNR].kappatype[i];
    if (!PisNULL(i)) {     
      if (type == STRSXP) {
	int endfor = NROW(i) * NCOL(i);
	for (int k=0; k<endfor; k++) FREE(PCHAR(i)[k]);
      } else if (isRObject(type)) {
	sexp_type *S = PSEXP(i);
	if (S->Delete) R_ReleaseObject(S->sexp);	
      } else if (type >= LISTOF) {
	LIST_DELETE((listoftype **) (cov->px + i));
      } // else nothing to free inside
      PFREE(i);
      cov->ncol[i] = cov->nrow[i] = SIZE_NOT_DETERMINED; // ==0
    }
  }

  MPPPROPERTIES_DELETE(&(cov->mpp));

  if (cov->ownkappanames != NULL) {
    int kappas = DefList[COVNR].kappas;
    for (j=0; j<kappas; j++) FREE(cov->ownkappanames[j]);
    UNCONDFREE(cov->ownkappanames);
  }

  QFREE;
 
  // important check in combination with above; can be easily removed or 
  // generalised  !!!!
 
  //MPPPROPERTIES_DELETE(&(cov->mpp));

  cov->prevloc = NULL;
  LOC_DELETE(&(cov->ownloc)); // OK ist es nicht, siehe oben

  if (cov->key != NULL) {
    //    TREE0(cov); /* ja nicht PMI, da dies auf geloeschtes zugreift */
    COV_DELETE(&(cov->key), save);
    //printf("done\n");
  }
  if (cov->rf != NULL && cov->origrf) UNCONDFREE(cov->rf);

  ce_DELETE(&(cov->Sce));
  localCE_DELETE(&(cov->SlocalCE));
  approxCE_DELETE(&(cov->SapproxCE));
  direct_DELETE(&(cov->Sdirect));
  hyper_DELETE(&(cov->Shyper));  
  nugget_DELETE(&(cov->Snugget));

  plus_DELETE(&(cov->Splus));
  sequ_DELETE(&(cov->Ssequ));
  //  SPECTRAL_DELETE(&(cov->Sspectral));
  tbm_DELETE(&(cov->Stbm));
  br_DELETE(&(cov->Sbr));
  get_DELETE(&(cov->Sget));
  pgs_DELETE(&(cov->Spgs));
  fctn_DELETE(&(cov->Sfctn));
  set_DELETE(&(cov->Sset));
  polygon_DELETE(&(cov->Spolygon));
  rect_DELETE(&(cov->Srect));
  dollar_DELETE(&(cov->Sdollar));
  gatter_DELETE(&(cov->Sgatter));
  earth_DELETE(&(cov->Searth));
  extra_DELETE(&(cov->Sextra));
  Ext_solve_DELETE(&(cov->Ssolve));
  biwm_DELETE(&(cov->Sbiwm));
  bistable_DELETE(&(cov->Sbistable));
  
  scatter_DELETE(&(cov->Sscatter));
  mcmc_DELETE(&(cov->Smcmc));
  //  SELECT_DELETE(&(cov->Sselect));
  likelihood_DELETE(&(cov->Slikelihood));  
  covariate_DELETE(&(cov->Scovariate));
  bubble_DELETE(&(cov->Sbubble));
  mle_DELETE(&(cov->Smle));
  model_DELETE(&(cov->Smodel), save);
  trend_DELETE(&(cov->Strend));
  gen_DELETE(&(cov->Sgen));
 
  simu_storage *simu = &(cov->simu);
  simu->active = simu->pair = false;
  simu->expected_number_simu = 0;
  if (cov->base!=NULL && cov->base->error_causing_cov == cov) {
    cov->base->error_causing_cov = save;
    if (save != NULL) STRCPY(save->err_msg, cov->err_msg);
  }

  UNCONDFREE(*Cov);
  *Cov = NULL;
}

void COV_DELETE_WITHOUT_LOC(model **Cov, model *save) { 
  model *cov = *Cov;
  int i,
    nsub = DefList[COVNR].maxsub;

  if (cov == NULL) {
    warning("*cov is NULL in COV_DELETE");
    return;
  }

  assert(cov != NULL);

  for (i=0; i<MAXPARAM; i++) { // cov->sub[i] // seit 1.10.07: luecken erlaubt
    //                         bei PSgen !
    if (cov->kappasub[i] != NULL) {
      //printf("del kappa %d ", i);
      COV_DELETE_WITHOUT_LOC(cov->kappasub + i, save);
    }
  }

  for (i=0; i<nsub; i++) { // cov->sub[i] // seit 1.10.07: luecken erlaubt
    //                         bei PSgen !
    // Achtung nach nsub koennen "blinde" Untermodelle kommen, die
    // nicht geloescht werden duerfen!!
    if (cov->sub[i] != NULL) {
      //printf("del sub %d ", i);
      COV_DELETE_WITHOUT_LOC(cov->sub + i, save);
    }
  }
  COV_DELETE_WITHOUTSUB(Cov, save);
}

// COV_DELETE(
void COV_DELETE_(model **Cov, model *save) { 
  if (*Cov == NULL) return; // 13.1.21 neu
  assert(*Cov != NULL); // TO DO 13.1.21 ehemals -> if ( != NULL) streichen
  model *cov = *Cov;
  //  PMI0(cov);  printf("COV_DELETE %.50s top=%d\n", NAME(cov), cov->calling == NULL);
  if (cov->calling == NULL) LOC_DELETE(&(cov->prevloc)); // OK
  COV_DELETE_WITHOUT_LOC(Cov, save);
}


void COV_ALWAYS_NULL(model *cov) {
  //  printf("\n\n\n\nCOVALWAYSNULL %s %d\n", NAME(cov), cov->base==NULL); 
  if (cov->base != NULL) cov->zaehler = ++(cov->base->zaehler);
  else cov->zaehler=0;
  
  cov->calling = NULL;
  cov->prevloc = cov->ownloc = NULL;
  cov->checked = false;

  cov->key = NULL;
  cov->rf = NULL;

  cov->mpp.mM = cov->mpp.mMplus = NULL;
  cov->mpp.moments = UNSET;
  cov->err_level = 0;
  cov->err = NOERROR;
  STRCPY(cov->err_msg, "<unknown>");
  
  cov->Sce = NULL;
  cov->SlocalCE = NULL;
  cov->SapproxCE = NULL;
  cov->Sdirect = NULL;
  cov->Shyper = NULL;
  cov->Snugget = NULL;
  cov->Splus = NULL;
  cov->Ssequ = NULL;
  cov->Stbm = NULL;
  cov->Sbr = NULL;
  cov->Sget = NULL;
  cov->Spgs = NULL;
  cov->Sfctn = NULL;
  cov->Sset = NULL;
  cov->Spolygon = NULL;
  cov->Srect = NULL;
  cov->Sdollar = NULL;
  cov->Sgatter = NULL;
  cov->Searth = NULL;
  cov->Sextra = NULL;
  cov->Ssolve = NULL;
  cov->Sbiwm = NULL;
  cov->Sbistable = NULL;
  cov->Sscatter = NULL;
  cov->Smcmc = NULL;
  //cov->Sselect = NULL;
  cov->Slikelihood = NULL;
  cov->Sbistable = NULL;
  cov->Scovariate = NULL;
  cov->Sbubble = NULL;
  cov->Smle = NULL;
  cov->Smodel = NULL;
  cov->Strend = NULL;
  cov->Sgen = NULL;

  cov->fieldreturn = falsch;
  cov->origrf = false;
  cov->initialised = false;
  cov->DallowedDone = false;
  cov->IallowedDone = false;
}


void SYSTEM_NULL(system_type *sys, int len) {
  for (int i=0; i<len; i++) {
    set_cumxmit(sys, i, UNSET);
    set_logdim(sys, i, UNSET);
    set_xdim_blank(sys, i, UNSET);
    set_maxdim(sys, i, UNSET);
    set_nri(sys, i, UNSET);
    set_last(sys, i, UNSET);
    set_type(sys, i, BadType);
    set_dom(sys, i, DOMAIN_MISMATCH);
    set_iso(sys, i, ISO_MISMATCH);    
  }
}


ALL_NULL(simu)
void COV_NULL(model *cov, KEY_type *base) {
  //printf("\n\n\n\nCCOVNULL %s %d\n", NAME(cov), base==NULL); 
  MEMSET(cov, 0, sizeof(model));
  if (base != NULL) cov->zaehler = ++(base->zaehler);// 0 never appears
  else cov->zaehler=0; // 0 is reserved also for use_external_boxcox
	 
  cov->mpp.moments = UNSET;

  set_nr(OWN, UNSET);
  cov->variant = UNSET;
  cov->user_given = ug_internal; // to know which models are given by the
  // user and which ones have been added internally

  cov->frame = BadType;
  cov->method = Forbidden;
  SYSTEM_NULL(PREV, MAXSYSTEMS);
  SYSTEM_NULL(GATTER, MAXSYSTEMS);
  SYSTEM_NULL(OWN, MAXSYSTEMS);
  VDIM0 = VDIM1 = UNSET;
  

  cov->ownkappanames = NULL;
  cov->logspeed = RF_NA;
  cov->full_derivs = cov->rese_derivs = UNSET;
  cov->ptwise_definite = pt_undefined;

  cov->monotone = MON_UNSET; 
  //  cov->total_n = UNSET;
  //  cov->total_sum = 0.0;
  // cov->diag = cov->semiseparatelast  = cov->separatelast = 
  cov->hess = NOT_IMPLEMENTED;
  int i;
  for (i=0; i<=Nothing; i++) cov->pref[i] = PREF_BEST;
  // mpp und extrG muessen auf NONE sein !
  for (; i<=Forbidden; i++) cov->pref[i] = PREF_NONE;

  MPPPROPERTIES_NULL(&(cov->mpp));
  simu_NULL(&(cov->simu));
  // cov->localcov=NULL;
}


ALL_NULL(FFT)
void FFT_destruct(FFT_storage *FFT)
{
  FREE(FFT->iwork);
  FREE(FFT->work);
  FFT_NULL(FFT);
}


ALL_NULL(ce)
void ce_DELETE(ce_storage **S) {
  ce_storage *x = *S;
  if (x != NULL) {
    int l,       
      vdim = x->vdim,      
      vdimSQ = vdim * vdim;
    if (x->c!=NULL) {
      for(l=0; l<vdimSQ; l++) FREE(x->c[l]);
      UNCONDFREE(x->c);
    }
    if (x->d!=NULL) {
      for(l=0; l<vdim; l++) FREE(x->d[l]);
      UNCONDFREE(x->d);
    }
#ifdef CE_PARALLEL
    for (int i=0; i<vdimSQ; i++) FFT_destruct(x->FFT + i);
    FREE(x->FFT);
#else
    FFT_destruct(&(x->FFT));
#endif
    
    FREE(x->aniso);
    FREE(x->gauss1);
    FREE(x->gauss2);
    UNCONDFREE(*S);
  }
}


ALL_NULL(localCE)
void localCE_DELETE(localCE_storage**S) 
{
  localCE_storage* x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}


ALL_NULL(approxCE)
void approxCE_DELETE(approxCE_storage **S) {
  approxCE_storage* x = * S;
  if (x != NULL) {
    FREE(x->idx);
    UNCONDFREE(*S);
  }
}


ALL_NULL(direct)
void direct_DELETE(direct_storage  ** S) {
  direct_storage *x = *S;
  if (x!=NULL) {
    FREE(x->G);
    UNCONDFREE(*S);
  }
}

void hyper_NULL(hyper_storage* x) { if (x == NULL) return; }
void hyper_DELETE(hyper_storage  **S) {
  hyper_storage *x = *S; 
  if (x != NULL) {
    UNCONDFREE(*S);
  }
}

ALL_NULL(nugget)
void nugget_DELETE(nugget_storage ** S)
{
  nugget_storage *x = *S;
  if (x != NULL) {
    FREE(x->pos);
    FREE(x->red_field);
    FREE(x->reduced_dim);
    FREE(x->datapos);
    FREE(x->prod_dim);
    FREE(x->index);
    UNCONDFREE(*S);
  }
}

ALL_NULL(plus)
void plus_DELETE(plus_storage ** S){
  plus_storage *x = *S;
  if (x != NULL) {
    UNCONDFREE(*S);
  }
}

ALL_NULL(sequ)
void sequ_DELETE(sequ_storage ** S){
  sequ_storage *x = *S;
  if (x!=NULL) {
    FREE(x->U11); 
    FREE(x->U22); 
    FREE(x->MuT); 
    FREE(x->G);
    FREE(x->Cov21);
    FREE(x->Inv22);
    FREE(x->res0);
    UNCONDFREE(*S);
  }
}


ALL_NULL(trend)
void trend_DELETE(trend_storage ** S) {
  trend_storage *x = *S;
  if (x!=NULL) {
    FREE(x->x);
    FREE(x->xi);
    FREE(x->evalplane);
    FREE(x->powmatrix);    
    UNCONDFREE(*S);
  }
}


void tbm_NULL(tbm_storage* x) {
  if (x == NULL) return;
  x->ce_dim =  0;
  x->err = ERRORFAILED;
}

void tbm_DELETE(tbm_storage **S) {
  tbm_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}


ALL_NULL(br)
void BRTREND_DELETE(double **BRtrend, int trendlen) {
  int j;
  if (BRtrend == NULL) return;
  for (j=0; j<trendlen; j++) FREE(BRtrend[j]);
}

void br_DELETE(br_storage **S) {
  br_storage *sBR = *S;  
  if (sBR != NULL) {
    assert(sBR->nr);
    
    int i;
    if (sBR->trend != NULL) {
      BRTREND_DELETE(sBR->trend, sBR->trendlen); 
      UNCONDFREE(sBR->trend);
    }

    if (sBR->nr == BRSHIFTED_INTERN || sBR->nr == BRSHIFTED_USER) {
      FREE(sBR->shift.loc2mem);
      FREE(sBR->shift.mem2location);
      FREE(sBR->shift.locindex);
    } else if (sBR->nr == BRMIXED_INTERN || sBR->nr == BRMIXED_USER) {
      if (sBR->m3.countvector != NULL) {
	for (i=0; i<sBR->m3.vertnumber; i++) FREE(sBR->m3.countvector[i]);
	UNCONDFREE(sBR->m3.countvector);
      }
      FREE(sBR->m3.loccentre);
      FREE(sBR->m3.logvertnumber);
      FREE(sBR->m3.lowerbounds);
      FREE(sBR->m3.suppmin);
      FREE(sBR->m3.suppmax);
      FREE(sBR->m3.areamatrix);
    } else if (sBR->nr == BRNORMED) {
      FREE(sBR->normed.current_prob);
      FREE(sBR->normed.current_cumprob);
      if (sBR->normed.C != NULL) {
	if (sBR->normed.do_not_delete_C) { FREE(sBR->normed.C[0]);
	} else {
	  int total = sBR->normed.total;
	  for (i=0; i<total; i++) FREE(sBR->normed.C[i]);
	}
	FREE(sBR->normed.C);
      }
      FREE(sBR->normed.dummyCi);    
    } else if (sBR->nr != BRORIGINAL_INTERN && sBR->nr != BRORIGINAL_USER) BUG;
         
    UNCONDFREE(*S);
  }
}


ALL_NULL_BUT(pgs, x->alpha = 1.0) // alpha!=1.0 only for optiz
void pgs_DELETE(pgs_storage **S) 
{
  pgs_storage *x = *S;
  if (x!=NULL) {
    // Huetchen
    FREE(x->v);
    FREE(x->y);
    FREE(x->xgr[0]);
    FREE(x->xgr);
   
    FREE(x->pos);
    FREE(x->min);
    FREE(x->max); 
    
    FREE(x->single);
    FREE(x->total);
    FREE(x->halfstepvector);    
    
    FREE(x->localmin);
    FREE(x->localmax);
    FREE(x->minmean);
    FREE(x->maxmean);
    
    // rf_interface.cc
    FREE(x->supportmin);
    FREE(x->supportmax);
    FREE(x->supportcentre);   
    FREE(x->own_grid_start);
    FREE(x->own_grid_step);
    FREE(x->own_grid_len);
    FREE(x->own_grid_cumsum);

    FREE(x->gridlen);
    FREE(x->end);
    FREE(x->start);
    FREE(x->delta);
    FREE(x->nx);

    FREE(x->xstart);    
    FREE(x->x);
    FREE(x->inc);    
 
    UNCONDFREE(*S);
  }
}



ALL_NULL(fctn) // alpha!=1.0 only for optiz
void fctn_DELETE(fctn_storage **S) 
{
  fctn_storage *x = *S;
  if (x!=NULL) {
    // grid
    FREE(x->start);
    FREE(x->end);
    FREE(x->nx);
    FREE(x->x);
    FREE(x->xstart);    
    FREE(x->inc);    
 
    FREE(x->startny);
    FREE(x->endy);
    FREE(x->ny);
    FREE(x->y);
    FREE(x->ystart);
    FREE(x->incy);    
 
    FREE(x->cross);
    FREE(x->cum);
     
    //
    FREE(x->C0x);
    FREE(x->C0y);
   
    UNCONDFREE(*S);
  }
}




ALL_NULL(set)
void set_DELETE(set_storage **S) {
  set_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}




void polygon_NULL(polygon_storage* x) {
  if (x == NULL) return;
  x->vprim = NULL;
  x->vdual = NULL;
  x->n_vdual = x->n_vertex = x->n_v = 0;

  polygon *P = x->P;
  if (P == NULL) BUG;
  P->e = NULL;
  P->v = NULL;
  P->n = 0;
}
void polygon_DELETE(polygon_storage **S) 
{
  polygon_storage *x = *S;
  if (x != NULL) {
    if (x->vdual != NULL) {
      int i;
      for (i=0; i<x->n_vdual; i++) FREE(x->vdual[i]);
      UNCONDFREE(x->vdual);
    }
    FREE(x->vprim);
    if (x->P != NULL) {
      freePolygon(x->P);
      UNCONDFREE(x->P);
    }
  }
  UNCONDFREE(*S);
}


ALL_NULL(rect)
void rect_DELETE(rect_storage **S){
  rect_storage *x = *S;
  if (x!=NULL) {
    FREE(x->value);
    FREE(x->weight);
    FREE(x->tmp_weight);
    FREE(x->right_endpoint);
    FREE(x->ysort);
    FREE(x->z);
    FREE(x->squeezed_dim);
    FREE(x->asSign);
    FREE(x->idx);
    UNCONDFREE(*S);
  }
}


ALL_NULL(dollar)
void dollar_DELETE(dollar_storage **S) 
{
  dollar_storage *x = *S;
  if (x!=NULL) {
    FREE(x->save_aniso);
    FREE(x->inv_aniso);
    FREE(x->len);
    FREE(x->sd);
    FREE(x->total);
    FREE(x->cumsum);
    FREE(x->proj);
    UNCONDFREE(*S);
  }
}

ALL_NULL(gatter)
void gatter_DELETE(gatter_storage **S) 
{
  gatter_storage *x = *S;
  if (x!=NULL) {
#ifdef DO_PARALLEL
    assert(x->z == NULL);
    assert(x->z1 == NULL);
#else    
    FREE(x->z);
    FREE(x->z1);
#endif    
    UNCONDFREE(*S);
  }
}


ALL_NULL(earth)
void earth_DELETE(earth_storage **S) 
{
  earth_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}


ALL_NULL(extra)
void extra_DELETE(extra_storage **S) 
{
  extra_storage *x = *S;
  if (x!=NULL) {
#ifdef DO_PARALLEL
    assert(x->a1 == NULL);
    assert(x->a2 == NULL);
    assert(x->a3 == NULL);
    assert(x->b1 == NULL);
    assert(x->b2 == NULL);
    assert(x->c1 == NULL);
    assert(x->c2 == NULL);
    assert(x->i1 == NULL);
    //   assert(x->i2 == NULL);
    assert(x->j1 == NULL);
    // assert(x->j2 == NULL);
    assert(x->k1 == NULL);
    assert(x->k2 == NULL);
#else
    FREE(x->a1);
    FREE(x->a2);
    FREE(x->a3);
    FREE(x->b1);
    FREE(x->b2);
    FREE(x->c1);
    FREE(x->c2);
    FREE(x->i1);
    //   FREE(x->i2);
    FREE(x->j1);
    // FREE(x->j2);
    FREE(x->k1);
    FREE(x->k2);
#endif    
    UNCONDFREE(*S);
  }
}


ALL_NULL(biwm)
void biwm_DELETE(biwm_storage **S) 
{
  biwm_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}


ALL_NULL(bistable)
void bistable_DELETE(bistable_storage **S)
{
  bistable_storage *x = *S;
  if (x!=NULL) {
    UNCONDFREE(*S);
  }
}


ALL_NULL(scatter)
void scatter_DELETE(scatter_storage **S) {
  scatter_storage *x = *S;
  if (x!=NULL) {
    FREE(x->min);
    FREE(x->max);
    FREE(x->step);
    UNCONDFREE(*S);
  }
}

ALL_NULL(mcmc)
void mcmc_DELETE(mcmc_storage **S) {
  mcmc_storage *x = *S;
  if (x!=NULL) {
    FREE(x->pos);
    FREE(x->deltapos);
    FREE(x->propdelta);
    FREE(x->proposed);
    UNCONDFREE(*S);
  }
}


ALL_NULL_BUT(get, x->param_nr = UNSET)
void get_DELETE(get_storage **S){
  get_storage *x = *S;
  if (x != NULL) {
    FREE(x->idx);
    UNCONDFREE(*S);
  }
}



void spectral_NULL(spectral_storage *x) {
  if (x == NULL) return;
  MEMSET(x, 0, sizeof(spectral_storage));
  x->phistep2d = x->phi2d = x->prop_factor = RF_NA;   
  x->nmetro = UNSET;
  x->sigma = (double) UNSET;
  for (int d=0; d<MAXTBMSPDIM; d++) {
    x->E[d] = x->sub_var_cum[d] = RF_NA;
  }
}

/*void spectral_DELETE(spectral_storage **S) { 
  spectral_storage *x = *S;
  if (x!=NULL) {
    // do NOT delete cov --- only pointer
    // spectral_storage *x; x =  *((spectral_storage**)S);
    UNCONDFREE(*S);   
  }
}
*/


void gen_NULL(gen_storage *x) {
  if (x == NULL) return; 
  // for (d=0; d<MAXMPPDIM; d++) 
  //   x->window.min[d] = x->window.max[d] = x->window.centre[d] = RF_NA;
  spectral_NULL(&(x->Sspectral));
  x->check = x->dosimulate = true;
}

void gen_DELETE(gen_storage **S) {
  gen_storage *x = *S;
  if (x != NULL) {
    //spectral_DELETE(&(x->Sspectral));
    UNCONDFREE(*S);
  }
}


void likelihood_facts_NULL(likelihood_facts *x) {
  if (x == NULL) return;
  x->varmodel = model_undefined;
  x->Var = NULL; // DO NOT FREE
  x->Matrix = x->pt_variance = NULL; // DO NOT FREE
  x->globalvariance = x->trans_inv = x->isotropic = false;
  x->newxdim = x->neffect = 0;  
  for (int i=0; i<MAX_LIN_COMP; i++) {
    x->effect[i] = 0;
    x->nas[i] = 0;
  }
}
void likelihood_facts_DELETE(likelihood_facts *x) {
  // nicht: Var, ptvariance !
  FREE(x->Matrix);
  //  x->Var = NULL; // DO NOT FREE
  //  x->Matrix = x->pt_variance = NULL; // DO NOT FREE
}


ALL_NULL_BUT(likelihood, x->facts.varmodel = model_undefined)
void likelihood_DELETE(likelihood_storage **S) {
  likelihood_storage *x = *S;
  if (x != NULL) {
    LIST_DELETE(&(x->datasets));
    if (x->X != NULL) for (int i=0; i<x->sets; i++) FREE(x->X[i]);
    FREE(x->X);
    if (x->YhatWithoutNA != NULL) for (int i=0; i<x->sets; i++) 
				    FREE(x->YhatWithoutNA [i]);
    FREE(x->YhatWithoutNA);
    FREE(x->XCY);
    FREE(x->XtX);
    FREE(x->XitXi);
    FREE(x->C);
    FREE(x->CinvXY);
    FREE(x->Cwork);
    FREE(x->Xwork);
    FREE(x->matrix);
    FREE(x->betavec);
    FREE(x->work);
    //FREE(x->Yhat);
    FREE(x->where_fixed);
    FREE(x->sumY);
    FREE(x->data_nas);  
    int end = x->cum_n_betas[x->fixedtrends];
    for (int i=0; i<end; i++) FREE(x->betanames[i]);    
    likelihood_facts_DELETE(&(x->facts)); 
    UNCONDFREE(*S);
  }
}

/*
  likelihood_storage * likelihood_CREATEXXX(int n) { /// obsolete
  likelihood_storage * x =
  (likelihood_st orage*) MAL LOC(sizeof(likelihood_storage));
  likelihood _NULL(x);
  x->datasets = LIST_CREATE(n, LISTOF + REALSXP);
  // x->withoutdet = LIST_CREATE(n, LISTOF + REALSXP);
  if ((x->X = (double**) CALLOC(n, sizeof(double*))) == NULL ||
  (x->XtX = (double*) CALLOC(n, sizeof(double))) == NULL ||
  (x->XitXi = (double*) CALLOC(n, sizeof(double))) == NULL ||
  (x->CinvXY = (double*) CALLOC(n, sizeof(double))) == NULL ||
  (x->data_nas = (bool*) CALLOC(n, sizeof(bool))) == NULL)
  err or("Memory allocation error for 'likelihood'");
  x->sets = n;
  }
*/

// #define PRNT printf("null: %lu\n", );

ALL_NULL_BUT(covariate, x->matrix_err = MATRIX_NOT_CHECK_YET )
void covariate_DELETE(covariate_storage **S) {
  covariate_storage *x = *S;
  if (x != NULL) {
    if (x->loc != NULL) LOC_DELETE(&(x->loc)); // OK da nur Parameter
    FREE(x->x);
    FREE(x->nrow);
    UNCONDFREE(*S);
   }
}


ALL_NULL(bubble)
void bubble_DELETE(bubble_storage **S) {
  bubble_storage *x = *S;
  if (x != NULL) {
    FREE(x->tau);
    FREE(x->rank);
    FREE(x->start);
    FREE(x->end);
    UNCONDFREE(*S);
  }
} 

ALL_NULL(mle)
void mle_DELETE(mle_storage **S) {
  mle_storage *x = *S;
  if (x != NULL) {
    UNCONDFREE(*S);
  }
}

ALL_NULL(model)
void model_DELETE(model_storage **S, model *save) {
  model_storage *x = *S;
  if (x != NULL) {
     COV_DELETE(&(x->vario), save);
     COV_DELETE(&(x->orig), save);
     COV_DELETE(&(x->get_cov), save);
     COV_DELETE(&(x->remote), save);
     for (int i=0; i<MAXSUB; i++) COV_DELETE(x->keys + i, save);

     if (x->cov != NULL) { // ehem. pgs_storage
       model *dummy = x->cov;
       if (x->cov->Smodel != NULL && x->cov->Smodel->cov != NULL && 
	   x->cov->Smodel->cov->Smodel == x) {
	 x->cov->Smodel->cov = NULL;// um "Unendlich"-Schleife zu vermeiden!
       } else {
	 assert(x->cov->Smodel == NULL || x->cov->Smodel->cov == NULL || 
		x->cov->Smodel->cov->Smodel == x);
       }
      x->cov = NULL; // Sicherheitshalber
      COV_DELETE(&(dummy), save);
     }     
    UNCONDFREE(*S);
  }
}

