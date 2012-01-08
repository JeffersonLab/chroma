// $Id: invcg2_timing_hacks_3.cc,v 3.1 2007-02-22 21:11:50 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "invcg2_timing_hacks_2.h"


// This is a hack version of invcg2 designed to allow timing.
// In particular stuff that would make the CG converge has been
// DISABLED... The length of the CG is completely controlled by 
// the INPUT n_count

//! Conjugate-Gradient (CGNE) algorithm for a generic Linear Operator
/*! \ingroup invert
 * This subroutine uses the Conjugate Gradient (CG) algorithm to find
 * the solution of the set of linear equations
 *
 *   	    Chi  =  A . Psi
 *
 * where       A = M^dag . M
 *
 * Algorithm:

 *  Psi[0]  :=  initial guess;    	       Linear interpolation (argument)
 *  r[0]    :=  Chi - M^dag . M . Psi[0] ;     Initial residual
 *  p[1]    :=  r[0] ;	       	       	       Initial direction
 *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;      Converged?
 *  FOR k FROM 1 TO MaxCG DO    	       CG iterations
 *      a[k] := |r[k-1]|**2 / <Mp[k],Mp[k]> ;
 *      Psi[k] += a[k] p[k] ;   	       New solution vector
 *      r[k] -= a[k] M^dag . M . p[k] ;        New residual
 *      IF |r[k]| <= RsdCG |Chi| THEN RETURN;  Converged?
 *      b[k+1] := |r[k]|**2 / |r[k-1]|**2 ;
 *      p[k+1] := r[k] + b[k+1] p[k];          New direction
 *
 * Arguments:
 *
 *  \param M       Linear Operator    	       (Read)
 *  \param chi     Source	               (Read)
 *  \param psi     Solution    	    	       (Modify)
 *  \param RsdCG   CG residual accuracy        (Read)
 *  \param MaxCG   Maximum CG iterations       (Read)
 *  \param n_count Number of CG iteration      (Write)
 *
 * Local Variables:
 *
 *  p   	       Direction vector
 *  r   	       Residual vector
 *  cp  	       | r[k] |**2
 *  c   	       | r[k-1] |**2
 *  k   	       CG iteration counter
 *  a   	       a[k]
 *  b   	       b[k+1]
 *  d   	       < p[k], A.p[k] >
 *  Mp  	       Temporary for  M.p
 *
 * Subroutines:
 *                             +               
 *  A       Apply matrix M or M  to vector
 *
 * Operations:
 *
 *  2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns )
 */

void InvCG2EvenOddPrecWilsLinOpTHack(const WilsonDslash &D,
			   const LFerm& chi,
			   LFerm& psi,
			   const LScal& mass,
			   const LScal& RsdCG, 
			   int MaxCG, 
			   int& n_count)
{
  // Always work on checkerboard 1. -- Chroma convention for even odd prec solver
  const Subset&  s = rb[1];

  // Length of vectors; -- subset length times Ns (Nc*Complex taken care of below)
  int n_3vec = (s.end() - s.start() + 1);

  // Useful Coefficient
  REAL mone = (REAL)(-1);

  // Stuff that usually gets packaged.
  // Use LScal-s for scalars to be able to use a*x + y notation
  LScal fact;
  AT_REAL(fact) = (REAL)(Nd + AT_REAL(mass));

  LScal invfact;
  AT_REAL(invfact)=(REAL)1/AT_REAL(fact);

  LScal mquarterinvfact;
  AT_REAL(mquarterinvfact) = -0.25*AT_REAL(invfact);

  //  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));

  // Initial norm -- could be passed in I suppose
  REAL chi_sq =(REAL)AT_REAL(norm2(chi,s));
  //QDPIO::cout << "chi_norm = " << sqrt(chi_sq) << endl;

  // ( Target residuum * || chi || )^2. Ignored in this test
  REAL rsd_sq = (AT_REAL(RsdCG) * AT_REAL(RsdCG)) * chi_sq;

  //                                            +
  //  r[0]  :=  Chi - A . Psi[0]    where  A = M  . M
    
  //                       +
  //  r  :=  [ Chi  -  M(u)  . M(u) . psi ]
  LFerm  r, mp, mmp;

  LFerm tmp1, tmp2, tmp3;

  // Apply M^{dag}M to psi
  
  // M -- cb issues taken care of by D.apply 
  // could save on function call overhead here by calling
  // Pete's routines directly

  D.apply(tmp1, psi, PLUS, 0);
  D.apply(tmp2, tmp1,PLUS , 1);

  // One day there will be a single routine for these
  tmp3[rb[1]] = mquarterinvfact*tmp2;
  mp[rb[1]] = fact*psi + tmp3;

  // Mdag
  D.apply(tmp1, mp, MINUS, 0);
  D.apply(tmp2, tmp1, MINUS, 1);
  tmp3[rb[1]] = mquarterinvfact*tmp2;
  mmp[rb[1]] = fact*psi + tmp3;

  // Do a vaxpy3 norm here to get r and LOCAL || r ||^2 together
  REAL cp;
  LDble sum; 
  //r[s] = chi - mmp; and || r[s] ||^2
  vaxpy3_norm( FIRST_ELEM(r,s), &mone, FIRST_ELEM(mmp,s), (REAL *)FIRST_ELEM(chi,s), 
	       n_3vec, &cp);

  // Do DP sum. Could save function call overhead here
  AT_REAL(sum) = (DOUBLE)cp; 
  QDPInternal::globalSum(sum);
  cp = (REAL)AT_REAL(sum); 

  //  p[1]  :=  r[0]
  LFerm p;
  p[s] = r;
  

  //QDPIO::cout << "InvCG: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq << endl;
  
  // Disable early termination
#if 0
  //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
  // cp is now a REAL so no need for the Chroma toBool() ism
  if ( cp  <=  rsd_sq )
  {
    n_count = 0;
    return;
  }
#endif

  //
  //  FOR k FROM 1 TO n_count  DO
  //
  REAL a, b, c, d;
  REAL *aptr;
  LScal as,bs; // For expressions
  REAL ma;

  // Go to n_count iters
  for(int k = 1; k <= n_count; ++k)
  {
    //  c  =  | r[k-1] |**2
    c = cp;

    //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
    //      	       	       	       	       	  +
    //  First compute  d  =  < p, A.p >  =  < p, M . M . p >  =  < M.p, M.p >
    //  Mp = M(u) * p
    D.apply(tmp1, p, PLUS, 0);
    D.apply(tmp2, tmp1,PLUS , 1);
    tmp3[rb[1]] = mquarterinvfact*tmp2;

    // replace this with a vaxpy norm
    // in the future will need a vaxpby_norm
    // mp[rb[1]] = fact*psi + tmp3;
    //  d = | mp | ** 2
    vaxpy3_norm(FIRST_ELEM(mp,s), 
		&AT_REAL(fact), 
		FIRST_ELEM(p,s), 
                FIRST_ELEM(tmp3,s), 
		n_3vec, 
		&d);
  
     // Global sum the local norm 
     AT_REAL(sum) = (DOUBLE)d;
     QDPInternal::globalSum(sum);
     d = (REAL)AT_REAL(sum);

    // Disable this. Set a = 1 to stop convergence
#if 0
    a = c/d;
#else 
    a = 1/1;
#endif

    //  Psi[k] += a[k] p[k]
    AT_REAL(as)  = a;   // For expressions below
    psi[s] += as * p;	/* 2 Nc Ns  flops */

    //  r[k] -= a[k] A . p[k] ;
    //      	       +            +
    //  r  =  r  -  M(u)  . Mp  =  M  . M . p  =  A . p
    D.apply(tmp1, mp, MINUS, 0);
    D.apply(tmp2, tmp1, MINUS, 1);
    tmp3[rb[1]] = mquarterinvfact*tmp2;
    mmp[rb[1]] = fact*mp + tmp3;

    // replace these with a vaxpy3 norm. 
    //r[s] -= a * mmp;
    //  cp  =  | r[k] |**2
    ma = -a;
    vaxpy3_norm(FIRST_ELEM(r,s), 
		&ma, 
		FIRST_ELEM(mmp,s), 
		FIRST_ELEM(r,s), n_3vec, &cp);

     // Global sum the local norm
     AT_REAL(sum)=(DOUBLE)cp;
     QDPInternal::globalSum(sum);
     cp = (REAL)AT_REAL(sum);


     //QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << endl;

    // Disable termination
#if 0
    // cp and rsd_sq are now REAL so no need for Chroma toBool() ism
    if ( cp  <=  rsd_sq )
    {
      n_count = k;
      return;
    }
#endif

    // Disable convergence
#if 0
    //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
    b = cp / c;
#else
    b = 1/1;
#endif

    // For expressions
    AT_REAL(bs) = b;
    //  p[k+1] := r[k] + b[k+1] p[k]
    p[s] = r + bs*p;	/* Nc Ns  flops */
  }
  // n_count = MaxCG;
  // QDP_error_exit("too many CG iterations: count = %d", n_count);
}

