
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2017 -- 2017 Martin Schlather

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



#ifndef RF_MACHINE_SCHLATHERS
#define RF_MACHINE_SCHLATHERS 1

#ifdef DO_PARALLEL
#include <omp.h>
#undef PRINTF
//#define PRINTF Rprintf
#define PRINTF if (omp_get_num_threads() > 1) { error("\n\nnever use Rprintf/PRINTF within parallel constructions!!\n\n"); } else Rprintf // OK
#endif

#ifdef SCHLATHERS_MACHINEX // SCHLATHERS_MACHINEX kein Tippfehler mit 
#ifndef RANDOMFIELDS_DEBUGGING
#undef assert
#define assert(X) if (X) {} else {PRINTF("\nassert failed: %s\n", #X); RFERROR("stopped.";}
#endif
#endif

//#undef assert
//#define assert(X) { BUG }

//#if defined ASSERT_GATTER
//#undef ASSERT_GATTER
//#endif
//
#undef ASSERT_GATTER
#define ASSERT_GATTER(Cov) assert(TrafoOK(Cov))


#define LOC_NOT_INITIALISED \
  SERR2("locations not initialised (%.50s line %d).", __FILE__, __LINE__)
#define LNRC(I,rc) __extension__ \
  ({assert(PLIST(I) != NULL); PLIST(I)->rc[SET_IDX(cov, I)];})



#define assert_sys(sys,X) {system_type *sys__sys##X = sys; assert(sys__sys##X != NULL);}
#define assert_cov(cov) { model *cov__cov = cov; assert(cov__cov != NULL); }
#define iCorrect(sys,X,s) {				\
    assert_sys(sys,X); assert(LASTi((sys)[0])>=0);	\
    assert(s>=0 && (s) <= LASTi((sys)[0]));		\
  }
#define PREVSYSOF(cov) __extension__({assert_cov(cov); (cov)->prev;})
#define GATTERSYSOF(cov)__extension__({assert_cov(cov); (cov)->gatter;})
#define SYSOF(cov) __extension__({assert_cov(cov); (cov)->own;})

//
#define LASTSYSTEM(sys) __extension__\
  ({assert_sys(sys,1); assert(LASTi((sys)[0])>=0); LASTi((sys)[0]); })


#define RESETLAST(sys) {assert_sys(sys,2); LASTi((sys)[0]) = 0; }
#define CUMXOHNE(sys,s) __extension__\
  ({iCorrect(sys,3,s); (s)<=0 ? 0 : CUMXMITi(sys[(s) - 1]); })
#define CUMXMIT(sys,s) __extension__\
  ({iCorrect(sys,4,s); CUMXMITi(sys[s]); })
#define TOTALXDIM(sys) __extension__({\
      assert_sys(sys,5); assert(CUMXMITi((sys)[LASTSYSTEM(sys)])>0);	\
      CUMXMITi(sys[LASTSYSTEM(sys)]); })
#define XDIM(sys,s) __extension__({assert_sys(sys,5a);			\
      assert((s)>=0 && (s) <= LASTSYSTEM(sys)); XDIMi((sys)[s]); })  // OK
#define LOGDIM(sys,s) __extension__\
  ({iCorrect(sys,6,s); LOGDIMi((sys)[s]); })
#define MAXDIM(sys,s) __extension__\
  ({iCorrect(sys,7,s); assert(MAXDIMi((sys)[s]) != 0); MAXDIMi((sys)[s]); })


#define ISO(sys,s) __extension__\
  ({iCorrect(sys,8,s); ISOi((sys)[s]); })
#define DOM(sys,s) __extension__\
  ({iCorrect(sys,9,s); DOMi((sys)[s]); })
#define SYSTYPE(sys,s) __extension__\
  ({iCorrect(sys,10,s); TYPEi((sys)[s]); })
#define ANYDIMOF(cov) __extension__({ \
       assert(OWNTOTALXDIM == PREVTOTALXDIM &&\
	     OWNTOTALXDIM == total_logicaldim(SYSOF(cov)) &&\
	     PREVLASTSYSTEM  == OWNLASTSYSTEM && \
	     Loctsdim(cov) == OWNTOTALXDIM); \
      OWNTOTALXDIM;})
//#define ANYDIMOF(cov)  __extension__({assert_cov(cov); assert(false); OWNTOTALXDIM; })
#define ANYOWNDIM __extension__\
  ({assert(OWNTOTALXDIM == OWNLOGDIM(OWNLASTSYSTEM)); OWNTOTALXDIM;})
#define SYSMODEL(sys) __extension__				\
  ({iCorrect(sys,11,0); NRi((sys)[0]); })
#define COPYSYS(to, from) __extension__\
  ({assert_sys(to,to3); assert_sys(from,from3); MEMCOPY(to, from, sizeof(system_type));})

#define  COPYALLSYSTEMS(to, from, keepnr) __extension__({		\
      assert_sys(to,to1); assert_sys(from,from1);				\
  system_type VARIABLE_IS_NOT_USED *to__=to, *from__=from;		\
  int nr_ = keepnr ? SYSMODEL(to) : MISMATCH;				\
  MEMCOPY(to, from, sizeof(Systems_type));				\
  if (keepnr) { set_nr(to, nr_); }					\
  /* muss hier mit SET_NR statt set_nr gearbeitet werden ?? */		\
   })

  
#define COPYALLSYSTEMS_COND(to, from, keepnr) __extension__({		\
      assert_sys(to,to2); assert_sys(from,from2);			\
    system_type *from__; from__ = to; from__ = from;			\
    bool gatter_set_ = LASTi(((cov)->gatter)[0]) >= 0;			\
    int nr_ = MISMATCH; if (gatter_set_ && keepnr) nr_= SYSMODEL(to);	\
    MEMCOPY(to, from__, sizeof(Systems_type));				\
    if (gatter_set_ && keepnr) { set_nr(to, nr_); }			\
    /* muss hier mit SET_NR statt set_nr gearbeitet werden ?? */	\
    })

#define SET_IDX(Cov, IDX) __extension__({\
      assert((Cov)!=NULL && (Cov)->nrow[IDX] >=0 && (Cov)->nrow[IDX]<=1000000);\
      ((Cov)->base->set % (Cov)->nrow[IDX]);})

#define CHECK(C,L,X,type,D,I,V,F) __extension__\
  ({assert((type)!=RandomType); check2X(C,L,X,type,D,I,V,F);})

#define CHECK_THROUGHOUT(C,P,type,D,I,V,F) __extension__		\
  ({assert((type)!=RandomType); check2Xthroughout(C,P,type,D,I,V,F);})

#define CHECK_NO_TRAFO(C,L,X,type,D,I,V,F) __extension__		\
  ({assert((type)!=RandomType); check2Xnotrafo(C,L,X,type,D,I,V,F);})



#define STRUCT(Cov, NM)  __extension__\
  ({ASSERT_GATTER(Cov); DefList[FIRSTGATTER].Struct(Cov, NM);})

// #define NN
#define PARAM(Cov, IDX) __extension__\
  ({assert(DefList[MODELNR(Cov)].kappatype[IDX] == REALSXP);  (Cov)->px[IDX];})
#define PARAMINT(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] == INTSXP);\
      ((int *) (Cov)->px[IDX]);})
#define PARAMCHAR(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] == STRSXP);	\
      (char **)(Cov)->px[IDX];})
//#define PARAMVEC(Cov, IDX) __extension__({assert(DefList[MODELNR(Cov)].kappatype[IDX] == VECSXP);  ((SEXP) (Cov)->px[IDX]);})
#define PARAMVEC(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] == VECSXP);  \
      ((sexp_type *) (Cov)->px[IDX])->sexp;})
#define PARAMENV(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] == ENVSXP);\
      ((sexp_type *) (Cov)->px[IDX]);})
#define PARAMLANG(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] == LANGSXP);\
      ((sexp_type *) (Cov)->px[IDX]);})
#define PARAMLIST(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] >= LISTOF);  \
      ((listoftype *) (Cov)->px[IDX]);})
#define LPARAM(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] == LISTOF + REALSXP); \
      assert(PARAMLIST(Cov, IDX) != NULL); \
      ((double *) (PARAMLIST(Cov, IDX)->lpx[SET_IDX(Cov, IDX)]));})
#define LPARAMINT(Cov, IDX) __extension__({\
      assert(DefList[MODELNR(Cov)].kappatype[IDX] == LISTOF + INTSXP); \
      assert(PARAMLIST(Cov, IDX) != NULL); \
      ((int *) (PARAMKLIST(Cov, IDX)->lpx[SET_IDX(Cov, IDX)]));})
//#define LISTLIST(Cov, IDX) __extension__({assert((Cov)->px != NULL); (freelistoftype *) (Cov)->px[IDX];})

#define PARAM0(Cov, IDX) __extension__\
  ({assert(PARAM(Cov, IDX) != NULL); PARAM(Cov, IDX)[0];})
#define PARAM0INT(Cov, IDX)  __extension__\
  ({assert(PARAMINT(Cov, IDX) != NULL); PARAMINT(Cov, IDX)[0];})
#define PARAM0CHAR(Cov, IDX)  __extension__\
  ({assert(PARAMCHAR(Cov, IDX) != NULL); PARAMCHAR(Cov, IDX)[0];})
#define LPARAM0(Cov, IDX) __extension__({\
      assert(LPARAM(Cov, IDX) != NULL); LPARAM(Cov, IDX)[0];})
#define LPARAM0INT(Cov, IDX)  __extension__({\
      assert(LPARAMINT(Cov, IDX) != NULL); LPARAMINT(Cov, IDX)[0];})


#define PCOPY(TO, FROM, IDX) 					\
  {if (!((TO) != NULL && (FROM) != NULL &&			\
	 (TO)->px[IDX] != NULL && (FROM)->px[IDX] &&IDX > 0 &&	\
	 (FROM)->ncol[IDX] == (TO)->ncol[IDX] &&		\
	 (FROM)->nrow[IDX] == (TO)->nrow[IDX] &&			\
	 (DefList[MODELNR(FROM)].kappatype[IDX] ==			\
	  DefList[MODELNR(TO)].kappatype[IDX]) &&			\
	 (DefList[MODELNR(FROM)].kappatype[IDX]==REALSXP ||		\
	  DefList[MODELNR(FROM)].kappatype[IDX]==INTSXP)  )) BUG;	\
    MEMCOPY((TO)->px[IDX], (FROM)->px[IDX],				\
	    ((FROM)->nrow[IDX]) * ((FROM)->ncol[IDX]) *			\
	    (DefList[MODELNR(FROM)].kappatype[IDX]==REALSXP ? sizeof(double) : \
	     DefList[MODELNR(FROM)].kappatype[IDX]==INTSXP ? sizeof(int) : \
	     -1));							\
  }

#define QcovALLOC(cov, nr) {\
    assert((cov)->q == NULL); 			\
    (cov)->qlen = (int) (nr);						\
    if (((cov)->q = (double*) CALLOC((int) (nr), sizeof(double))) == NULL) \
      RFERROR("memory allocation error for local memory");		\
  }

#define B2D(X) ((X) != true) && ((X) != false) ? ({BUG; 0.0;}) : (double) (X)
#define B2I(X) (int) (X)
#define E2D(X) (X)
#define E2I(X) (int) (X)
#define E2B(X) (bool) (X)
#define I2D(X) (double) (X)


#define RootOf(cov) __extension__({				\
      assert(cov != NULL && cov->root != NULL); \
      cov->root;})
#define KEYtypeOf(cov) __extension__({				\
      assert(cov != NULL && cov->base != NULL); \
      cov->base;})
#define KT_NAOK_RANGE  __extension__({ assert(KT != NULL); KT->naok_range;})
#define setKT_NAOK_RANGE(KT,X) { assert(KT != NULL); KT->naok_range = X;}
#define NAOK_RANGE KEYtypeOf(cov)->naok_range
#define set_NAOK_RANGE(X) {				\
    assert(cov != NULL && cov->base != NULL);		\
    KEYtypeOf(cov)->naok_range = (X);}


#undef LocP
#undef Loc
#undef Loctsdim
#undef Loctsdim
#undef LocxdimOZ
#undef Locspatialdim
#undef Loctotalpoints

#define Loc(cov) __extension__({			\
      assert(cov != NULL);					\
      assert((cov)->ownloc!=NULL || (cov)->prevloc != NULL);	\
      ((cov)->ownloc!=NULL ? LocLoc((cov)->ownloc, cov) :	\
       LocLoc((cov)->prevloc, cov));})
#define LocP(cov) __extension__({				\
      assert(cov != NULL);					\
      assert(((cov)->ownloc != NULL ? (cov)->ownloc : (cov)->prevloc)!=NULL); \
      ((cov)->ownloc != NULL ? (cov)->ownloc : (cov)->prevloc) ;	\
     })
#define Loctsdim(cov) __extension__({assert( (LocP(cov) != NULL ? Loc(cov)->timespacedim : 0) > 0); (LocP(cov) != NULL ? Loc(cov)->timespacedim : 0);})
#define LocxdimOZ(cov) __extension__({assert( (LocP(cov) != NULL ? Loc(cov)->xdimOZ : 0)); (LocP(cov) != NULL ? Loc(cov)->xdimOZ : 0);})
#define Locspatialdim(cov) __extension__({assert((LocP(cov) != NULL ? Loc(cov)->spatialdim : 0));(LocP(cov) != NULL ? Loc(cov)->spatialdim : 0);})
#define Loctotalpoints(cov) __extension__({assert((LocP(cov) != NULL ? Loc(cov)->totalpoints : 0)); (LocP(cov) != NULL ? Loc(cov)->totalpoints : 0);})


#endif
