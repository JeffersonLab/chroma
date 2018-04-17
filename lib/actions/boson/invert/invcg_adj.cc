/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "invcg_adj.h"

using namespace QDP::Hints;
#undef PAT
#ifdef PAT
#include <pat_api.h>
#endif

namespace Chroma 
{

  Double norm2_adj(const LatticeColorMatrix& C,const Subset& s){
    return Double(0.5)*real(sum(trace(adj(C)*C),s)) ;
  }

  DComplex inner_prod_adj(const LatticeColorMatrix& C,
			const LatticeColorMatrix& X,
			const Subset& s){
    return Double(0.5)*sum(trace(adj(C)*X),s) ;
  }
  

  
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
   *      Psi[k] += a[k] p[k] ;   	       New solution std::vector
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
   *  \return res    System solver results
   *
   * Local Variables:
   *
   *  p   	       Direction std::vector
   *  r   	       Residual std::vector
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
   *  A       Apply matrix M or M  to std::vector
   *
   * Operations:
   *
   *  2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns )
   */

  SystemSolverResults_t 
  InvCG_adj(const LinearOperator<LatticeColorMatrix>& M,
	      const LatticeColorMatrix& chi,
	      LatticeColorMatrix& psi,
	      const Real& RsdCG, 
	      int MaxCG)
  {
    START_CODE();

    const Subset& s = M.subset();

    SystemSolverResults_t  res;
    LatticeColorMatrix mp;                moveToFastMemoryHint(mp);
    LatticeColorMatrix p;                 moveToFastMemoryHint(p);
    moveToFastMemoryHint(psi,true);
    LatticeColorMatrix r;                 moveToFastMemoryHint(r);
    LatticeColorMatrix chi_internal;      moveToFastMemoryHint(chi_internal);

    chi_internal[s] = chi;

    QDPIO::cout << "InvCG_adj: starting" << std::endl;
    FlopCounter flopcount;
    flopcount.reset();
    StopWatch swatch;
    swatch.reset();
    swatch.start();

//  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));
    // implements the norm
    Double chi_sq = norm2_adj(chi_internal,s);
    flopcount.addSiteFlops(4*Nc*Nc,s);

#if 1
    QDPIO::cout << "chi_norm = " << sqrt(chi_sq) << std::endl;
#endif

    Double rsd_sq = (RsdCG * RsdCG) * chi_sq;

    //                                
    //  r[0]  :=  Chi - A . Psi[0]    
    
    M(mp, psi,PLUS);
    flopcount.addFlops(M.nFlops());

    r[s] = chi_internal - mp;
    flopcount.addSiteFlops(2*Nc*Nc,s);
    //  Cp = |r[0]|^2
    Double cp = norm2_adj(r, s);   	       	   /* 2 Nc Ns  flops */
    flopcount.addSiteFlops(4*Nc*Ns, s);


#if 1
    QDPIO::cout << "InvCG: k = 0  || r ||= " <<sqrt(cp) << std::endl;
#endif

    //  p[1]  :=  r[0]
    p[s] = r;
  



    //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
    if ( toBool(cp  <=  rsd_sq) )
    {
      res.n_count = 0;
      res.resid   = sqrt(cp);
      swatch.stop();
      flopcount.report("invcg_adj", swatch.getTimeInSeconds());
      revertFromFastMemoryHint(psi,true);
      END_CODE();
      return res;
    }

    //
    //  FOR k FROM 1 TO MaxCG DO
    //

    Double a, b, c, d;
  
    for(int k = 1; k <= MaxCG; ++k)
    {
#if 1
      QDPIO::cout << "InvCG: k = "<<k<<"  || r ||= " <<sqrt(cp) << std::endl;
#endif

      //  c  =  | r[k-1] |**2
      c = cp;

      //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
      //      	       	       	       	       	  +
      //  First compute  d  =  < p, A.p >  =  < p, M . M . p >  =  < M.p, M.p >
      //  Mp = M(u) * p
      M(mp, p,PLUS);  flopcount.addFlops(M.nFlops());

      //  d = p' m p
      // m is hermitian hence d is real
      d = real(inner_prod_adj(p,mp, s));  flopcount.addSiteFlops(4*Nc*Nc,s);

      //  r[k] -= a[k] A . p[k] ;
      //      	       +            +
      //  r  =  r  -  M(u)  . Mp  =  M  . M . p  =  A . p
 
      a = c/d;

      Real ar = a;

      r[s] -= ar * mp;
      flopcount.addSiteFlops(4*Nc*Nc, s);

      //  cp  =  | r[k] |**2
      cp = norm2_adj(r, s);    flopcount.addSiteFlops(4*Nc*Nc,s);

      //  Psi[k] += a[k] p[k]
      psi[s] += ar * p;    flopcount.addSiteFlops(4*Nc*Nc,s);



      //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;


//    QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << std::endl;

      if ( toBool(cp  <=  rsd_sq) )
      {
	res.n_count = k;
	res.resid   = sqrt(cp);
	swatch.stop();
	//	QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << std::endl;
	flopcount.report("invcg2", swatch.getTimeInSeconds());
	revertFromFastMemoryHint(psi,true);

	// Compute the actual residual
	{
	  M(mp, psi,PLUS);
	  Double actual_res = norm2_adj(chi - mp,s);
	  res.resid = sqrt(actual_res);
	}

	END_CODE();
	return res;
      }

      //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
      b = cp / c;
      Real br = b;

      //  p[k+1] := r[k] + b[k+1] p[k]
      p[s] = r + br*p;    flopcount.addSiteFlops(4*Nc*Nc,s);
    }
    res.n_count = MaxCG;
    res.resid   = sqrt(cp);
    swatch.stop();
    QDPIO::cerr << "Nonconvergence Warning" << std::endl;
    flopcount.report("invcg_adj", swatch.getTimeInSeconds());
    revertFromFastMemoryHint(psi,true);
    QDPIO::cerr << "too many CG iterations: count =" << res.n_count <<" rsd^2= " << cp << std::endl <<std::flush;

    END_CODE();
    return res;
  }


}  // end namespace Chroma
