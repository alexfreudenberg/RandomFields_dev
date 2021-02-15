

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

/// sysconf (_SC_NPROCESSORS_ONLN) // number of cores available
// int get_nprocs (void) // dito

#ifndef RFdef_H
#define RFdef_H 1


#include "intrinsics.h"
#include "basic.h"

//
////  1
//
//
// 1



#include <Basic_utils.h>
#define DO_PARALLEL_ALREADY_CONSIDERED 1 // verhindert dass in  Basic_utils.h
//                                          entsprechend gesetzt wird
#ifdef _OPENMP
  #ifdef SCHLATHERS_MACHINE
  #else
    #define DO_PARALLEL 1
    #define CE_PARALLEL 1 // for circulant embedding
//#define DO_TALLOC_SAVE 1
// brownresnick, tbm
  #endif
#else
   #ifdef DO_PARALLEL
     #undef DO_PARALLEL
   #endif
//   #ifdef DO_TALLOC_SAVE 
//    #undef DO_TALLOC_SAVE 
//  #endif
#endif

#define DO_TALLOC_SAVE 1

#ifdef DO_PARALLEL
//#undef DO_PARALLEL
#endif

//  assert((cov) != NULL && (cov)->base != NULL);			

#define SAVE_ERR(cov, N)	\
  (cov)->err = N;  \
  (cov)->base->error_causing_cov =				\
    (N) == NOERROR ? NULL					\
    : (cov)->base->error_causing_cov == NULL ? (cov)		\
    : (cov)->base->error_causing_cov; 

#define RETURN_ERR(err) { SAVE_ERR(cov, err); return err; }
#define RETURN_ERR_COV(cov, err) { SAVE_ERR(cov, err); return err; }

#define RETURN_NOERROR { cov->err = NOERROR;\
    cov->base->error_causing_cov = NULL; \
    return NOERROR; }

#define LOCAL_ERROR(N) SAVE_ERR(cov, N)
#define WHICH_ERRORSTRING cov->err_msg // wird benoetigt
#define ERROR_LOC ""

#define LOCAL_ERRLOC_MSG char ERRMSG[LENERRMSG];
#define LOCAL_ERRORSTRING




#endif
