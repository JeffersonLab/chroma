// $Id: bicgstab.h,v 1.2 2009-10-21 21:03:39 kostas Exp $
/*------------------------------------------------------------------------
  BICGSTAB algorithm for solving A*x=b with A a sparse non-symmetric matrix 
  Authors     : Abdou M. Abdel-Rehim, Kostas Orginos, Andreas Stathopoulos
                andreas@cs.wm.edu
  Last Updated: August, 28th, 2009.
--------------------------------------------------------------------------*/

#ifndef MY_BICGSTAB_H_
#define MY_BICGSTAB_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "qdp-lapack_Complex.h"
#include "qdp-lapack_numerical.h"
#include "qdp-lapack_numerical_private.h"

void MyBICGSTAB_C(int n, int lde, Complex_C *x, Complex_C *b, float tol, int maxiter, int *iter, float *res,
                   void (*matvec) (void *, void *, void *), void *params, float AnormEst, int ConvTestOpt, 
                   Complex_C *work, int *flag);


/* BICGSTAB to solve the linear system Ax=b 
   n (IN)         : the active size of the matrix A.
   lde(IN)        : Leading dimension of A.
   x  (IN/OUT)    : on input this std::vector has the inital guess and on output it has the solution.
   b  (IN)        : right-hand side.
   tol(IN)        : tolerance for convergence. At convergence, norm(r) <= tol*norm(b).
   maxiter(IN)    : allowed maximum number of iterations.
   iter (OUT)     : number of iterations till convergence.
   res (OUT)      : std::vector of length iter of residual norms at every iteration.
   matvec (IN)    : matrix-std::vector product function.
   params (IN)    : parameters needed by the matrix-std::vector function.
   AnormEst(IN)   : Estimate of the norm of A.Irrelevant if ConvTestOpt ==1.
   ConvTestOpt (IN): choice of how to determine convergence:
                       1 means norm(res) < tol*norm(b) 
                       2 means norm(res) < max(tol*norm(b), Machine_precision*(AnormEst*norm(x)+norm(b)))
                         with Machine_precision=1e-7
   work(IN)       : work array of dimension >= 6*lde.
   flag(OUT)      : exit status
                    0 means the system reached convergence on a number of iterations less than or equal maxiter.
                    1 means that the right hand side is zero and no iterations were performed and zero solution is returned.
                    2 one of bicgstab scalar quantities were zero.
                    3 no convergence after maximum number of iterations performed.

--------------------------------------------------------------------------------------------*/
void MyBICGSTAB_Z(int n, int lde, Complex_Z *x, Complex_Z *b, double tol, int maxiter, int *iter, double *res,
                void (*matvec) (void *, void *, void *), void *params, double AnormEst, int ConvTestOpt,
                Complex_Z *work, int *flag);

                 //same as MyBICGSTAB_C but with all double precision and now Machine_precision=1e-16;

#endif
