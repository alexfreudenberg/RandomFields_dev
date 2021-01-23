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

#include "errors_messages.h"

#ifndef rfutils_solve_H
#define rfutils_solve_H 1


typedef enum InversionMethod { 
  Cholesky, // 0
  SVD,  // 1
  Eigen, // 2
  Sparse, // 3
  NoInversionMethod, // 4, last user available method
  QR, // 5
  LU, // 6 currently not propagated
  NoFurtherInversionMethod, // 7, local values
  GPUcholesky,		    // 8
  direct_formula,           // 9
  Diagonal // 10, always last one!
} InversionMethod;
#define nr_InversionMethods ((int) Diagonal + 1)
#define nr_user_InversionMethods ((int) NoFurtherInversionMethod + 1)


#define SOLVE_METHODS 3
typedef // benoetigt
struct solve_storage {
  errorstring_type err_msg;
  InversionMethod method, newMethods[SOLVE_METHODS];
  usr_bool sparse;
  int
  //   SICH_n, MM_n, VT_n, U_n, D_n, 
  //    work_n, w2_n, lnz_n, w3_n, result_n,  nsuper, nnzlindx,
    size, actual_size, actual_pivot,
    nsuper,
    main_n, rhs_n, w2_n, U_n, D_n, w3_n, lnz_n, result_n, 
    
    *pivot_idx, pivot_idx_n, 
    *iwork, iwork_n, //eigen, svd, LU, spam
    *pivotsparse, pivotsparse_n, *xlnz,xlnz_n, //spam
    *snode,snode_n, *xsuper, xsuper_n,*invp,invp_n,   // spam
    *cols,cols_n, *rows,rows_n, *lindx, lindx_n, // spam
    *xja,  xja_n; // chol, eigen, spam
  double 
  *main, *rhs,// diagonal, general -- FORBIDDEN for further use
    *w2, // eigen, svd, LU, QR, pivot
    *U, // eigen, svd, pivot
    *D, // eigen, svd, cholesky, spam, pivot
    *w3, // spam, QR, svd
    *lnz, // spam
    *result,  // sqrtPosDefFree
    *to_be_deleted; 
} solve_storage;

#define SOLVE 0
#define MATRIXSQRT 1
#define DETERMINANT 2

#endif
