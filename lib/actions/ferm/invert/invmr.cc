// $Id: invmr.cc,v 3.0 2006-04-03 04:58:49 edwards Exp $

/*! \file
 *  \brief Minimal-Residual (MR) for a generic fermion Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invmr.h"

namespace Chroma {

/*! \ingroup invert
 * This subroutine uses the Minimal Residual (MR) algorithm to determine
 * the solution of the set of linear equations. Here we allow A to be nonhermitian.

 *   	    Chi  =  A . Psi 

 * Algorithm:

 *  Psi[0]                                      Argument
 *  r[0]    :=  Chi  -  A . Psi[0] ;            Initial residual
 *  IF |r[0]| <= RsdMR |Chi| THEN RETURN;       Converged?
 *  FOR k FROM 1 TO MaxCG DO                    MR iterations
 *      a[k-1]  := <A.r[k-1],r[k-1]> / <A.r[k-1],A.r[k-1]> ;
 *      ap[k-1] := MRovpar * a[k] ;             Overrelaxtion step
 *      Psi[k]  += ap[k-1] r[k-1] ;   	        New solution vector
 *      r[k]    -= ap[k-1] A . r[k-1] ;         New residual
 *      IF |r[k]| <= RsdMR |Chi| THEN RETURN;   Converged?

 * Arguments:

 *  \param M       Linear Operator    	       (Read)
 *  \param chi     Source	               (Read)
 *  \param psi     Solution    	    	       (Modify)
 *  \param RsdMR   MR residual accuracy        (Read)
 *  \param MRovpar Overrelaxation parameter    (Read)
 *  \param MaxCG   Maximum CG iterations       (Read)
 *  \param n_count Number of CG iteration      (Write)

 * Local Variables:

 *  r   	Residual vector
 *  cp  	| r[k] |**2
 *  c   	| r[k-1] |**2
 *  k   	MR iteration counter
 *  a   	a[k]
 *  d   	< M.r[k], M.r[k] >
 *  R_Aux       Temporary for  A.Psi
 *  Mr 	        Temporary for  A.r

 * Global Variables:

 *  MaxCG       Maximum number of MR iterations allowed
 *  RsdMR       Maximum acceptable MR residual (relative to source)
 *
 * Subroutines:
 *
 *  A           Apply matrix to vector
 *
 * Operations:
 *
 *  ???? 
 */

void InvCG2(const LinearOperator& A,
	    const LatticeFermion& chi,
	    LatticeFermion& psi,
	    const Real& RsdMR, 
	    const Real& MRovpar;
	    int MaxCG, 
	    int& n_count)
{
  START_CODE();

  const OrderedSubset& s = M.subset();

  LatticeFermion Ar;
  Complex a;
  DComplex c;
  Double d;
  int k;

  Real rsd_sq = (RsdMR * RsdMR) * Real(norm2(chi,s));
        
  /*  r[0]  :=  Chi - A . Psi[0] */
  /*  r  :=  A . Psi  */
  LatticeFermion r;
  r[s] = chi - A(psi, PLUS);
  
  /*  Cp = |r[0]|^2 */
  Double cp = norm2(r, s);                 /* 2 Nc Ns  flops */

  /*  IF |r[0]| <= RsdMR |Chi| THEN RETURN; */
  /*#erite(nml_out,*) "&Minimal_Residual_iterations"; */
  /*#printf("InvMRm: k = 0  cp = %g\n", cp); */
  /*#write(nml_out,*) "ClovInv: k = 0  cp = ", cp; */
  
  if ( toBool(cp  <=  rsd_sq) )
  {
    n_count = 0;
    END_CODE();
    return;
  }

  /*  FOR k FROM 1 TO MaxCG DO */
  k = 0;
  /*print *, "MaxCG = ", MaxCG; */
  /*print *, "(RsdMR * chi_norm)**2 = ", rsq_sq; */
  while( (k < MaxCG) && (toBool(cp > rsd_sq)) )
  {
    ++k;

    /*  a[k-1] := < A.r[k-1], r[k-1] >/ < A.r[k-1], A.r[k-1] > ; */
    /*  Ar = A * r  */
    Ar[s] = A(r, PLUS);

    /*  c = < A.r, r > */
    c = innerProduct(Ar, r, s);
    
    /*  d = | A.r | ** 2  */
    d = norm2(Ar, s);	                /* 2 Nc Ns  flops */

    /*  a = c / d */
    a = c / d;
    
    /*  a[k-1] *= MRovpar ; */
    a = a * MRovpar;

    /*  Psi[k] += a[k-1] r[k-1] ; */
    psi[s] += r * a;	                /* 2 Nc Ns  flops */

    /*  r[k] -= a[k-1] M . r[k-1] ; */
    r[s] -= Ar * a;                     /* 2 Nc Ns  flops */

    /*  cp  =  | r[k] |**2 */
    cp = norm2(r, s);                   /* 2 Nc Ns  flops */

    /*#  printf("InvMRm: k = %d  cp = %g\n", k, cp); */
    /*#  write(nml_out,*) "InvMRm: k = ",k,"  cp = ", cp; */
  }
  n_count = k;
  if ( n_count == MaxCG )
    QDP_error_exit("too many MR iterations: iters=%d", k);

  END_CODE();
  return;
}

}  // end namespace Chroma
