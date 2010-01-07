// $Id: invcg2.cc,v 3.7 2009-06-02 15:56:39 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invcg2.h"

using namespace QDP::Hints;
#undef PAT
#ifdef PAT
#include <pat_api.h>
#endif

namespace Chroma 
{

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
   *  \return res    System solver results
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

  template<typename T, typename RT>
  SystemSolverResults_t 
  InvCG2_a(const LinearOperator<T>& M,
	   const T& chi,
	   T& psi,
	   const RT& RsdCG, 
	   int MaxCG)
  {
    START_CODE();

    const Subset& s = M.subset();

    SystemSolverResults_t  res;
    T mp;                moveToFastMemoryHint(mp);
    T mmp;               moveToFastMemoryHint(mmp);
    T p;                 moveToFastMemoryHint(p);
    moveToFastMemoryHint(psi,true);
    T r;                 moveToFastMemoryHint(r);
    T chi_internal;      moveToFastMemoryHint(chi_internal);

    chi_internal[s] = chi;

    QDPIO::cout << "InvCG2: starting" << endl;
    FlopCounter flopcount;
    flopcount.reset();
    StopWatch swatch;
    swatch.reset();
    swatch.start();

//  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));
    Double chi_sq =  norm2(chi_internal,s);
    flopcount.addSiteFlops(4*Nc*Ns,s);

#if 0
    QDPIO::cout << "chi_norm = " << sqrt(chi_sq) << endl;
#endif

    Double rsd_sq = (RsdCG * RsdCG) * chi_sq;

    //                                            +
    //  r[0]  :=  Chi - A . Psi[0]    where  A = M  . M
    
    //                      +
    //  r  :=  [ Chi  -  M(u)  . M(u) . psi ]
    M(mp, psi, PLUS);
    M(mmp, mp, MINUS);
    flopcount.addFlops(2*M.nFlops());

    r[s] = chi_internal - mmp;
    flopcount.addSiteFlops(2*Nc*Ns,s);
    //  Cp = |r[0]|^2
    Double cp = norm2(r, s);   	       	   /* 2 Nc Ns  flops */
    flopcount.addSiteFlops(4*Nc*Ns, s);


#if 0
    QDPIO::cout << "InvCG: k = 0  || r ||= " <<sqrt(cp) << endl;
#endif

    //  p[1]  :=  r[0]
    p[s] = r;
  



    //  IF |r[0]| <= RsdCG |Chi| THEN RETURN;
    if ( toBool(cp  <=  rsd_sq) )
    {
      res.n_count = 0;
      res.resid   = sqrt(cp);
      swatch.stop();
      flopcount.report("invcg2", swatch.getTimeInSeconds());
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
      //  c  =  | r[k-1] |**2
      c = cp;

      //  a[k] := | r[k-1] |**2 / < p[k], Ap[k] > ;
      //      	       	       	       	       	  +
      //  First compute  d  =  < p, A.p >  =  < p, M . M . p >  =  < M.p, M.p >
      //  Mp = M(u) * p
      M(mp, p, PLUS);  flopcount.addFlops(M.nFlops());

      //  d = | mp | ** 2
      d = norm2(mp, s);  flopcount.addSiteFlops(4*Nc*Ns,s);

      //  r[k] -= a[k] A . p[k] ;
      //      	       +            +
      //  r  =  r  -  M(u)  . Mp  =  M  . M . p  =  A . p
      M(mmp, mp, MINUS);
      flopcount.addFlops(M.nFlops());

 
      a = c/d;

      RT ar = a;

      r[s] -= ar * mmp;
      flopcount.addSiteFlops(4*Nc*Ns, s);

      //  cp  =  | r[k] |**2
      cp = norm2(r, s);    flopcount.addSiteFlops(4*Nc*Ns,s);

      //  Psi[k] += a[k] p[k]
      psi[s] += ar * p;    flopcount.addSiteFlops(4*Nc*Ns,s);



      //  IF |r[k]| <= RsdCG |Chi| THEN RETURN;


//    QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << endl;

      if ( toBool(cp  <=  rsd_sq) )
      {
	res.n_count = k;
	res.resid   = sqrt(cp);
	swatch.stop();
	//	QDPIO::cout << "InvCG: k = " << k << "  cp = " << cp << endl;
	flopcount.report("invcg2", swatch.getTimeInSeconds());
	revertFromFastMemoryHint(psi,true);

	// Compute the actual residual
	{
	  M(mp, psi, PLUS);
	  M(mmp, mp, MINUS);
	  Double actual_res = norm2(chi - mmp,s);
	  res.resid = sqrt(actual_res);
	}

	END_CODE();
	return res;
      }

      //  b[k+1] := |r[k]|**2 / |r[k-1]|**2
      b = cp / c;
      RT br = b;

      //  p[k+1] := r[k] + b[k+1] p[k]
      p[s] = r + br*p;    flopcount.addSiteFlops(4*Nc*Ns,s);
    }
    res.n_count = MaxCG;
    res.resid   = sqrt(cp);
    swatch.stop();
    QDPIO::cerr << "Nonconvergence Warning" << endl;
    flopcount.report("invcg2", swatch.getTimeInSeconds());
    revertFromFastMemoryHint(psi,true);
    QDPIO::cerr << "too many CG iterations: count =" << res.n_count <<" rsd^2= " << cp << endl <<flush;

    END_CODE();
    return res;
  }


  //
  // Explicit versions
  //
  // Single precision
  SystemSolverResults_t 
  InvCG2(const LinearOperator<LatticeFermionF>& M,
	 const LatticeFermionF& chi,
	 LatticeFermionF& psi,
	 const Real& RsdCG, 
	 int MaxCG)
  {
#ifdef PAT
    int ierr = PAT_region_begin(20, "InvCG2Single");
#endif
    return InvCG2_a<LatticeFermionF,RealF>(M, chi, psi, RsdCG, MaxCG);
#ifdef PAT
    ierr = PAT_region_end(20);
#endif

  }

  // Double precision
  SystemSolverResults_t 
  InvCG2(const LinearOperator<LatticeFermionD>& M,
	 const LatticeFermionD& chi,
	 LatticeFermionD& psi,
	 const Real& RsdCG, 
	 int MaxCG)
  {
#ifdef PAT
    int ierr=PAT_region_begin(21, "InvCG2Double");
#endif
    return InvCG2_a<LatticeFermionD, RealD>(M, chi, psi, RsdCG, MaxCG);
#ifdef PAT
    ierr= PAT_region_end(21);
#endif
  }

  // Single precision
  SystemSolverResults_t 
  InvCG2(const LinearOperator<LatticeStaggeredFermionF>& M,
	 const LatticeStaggeredFermionF& chi,
	 LatticeStaggeredFermionF& psi,
	 const Real& RsdCG, 
	 int MaxCG)
  {
    return InvCG2_a<LatticeStaggeredFermionF,RealF>(M, chi, psi, RsdCG, MaxCG);
  }

  // Double precision
  SystemSolverResults_t 
  InvCG2(const LinearOperator<LatticeStaggeredFermionD>& M,
	 const LatticeStaggeredFermionD& chi,
	 LatticeStaggeredFermionD& psi,
	 const Real& RsdCG, 
	 int MaxCG)
  {
    return InvCG2_a<LatticeStaggeredFermionD,RealD>(M, chi, psi, RsdCG, MaxCG);
  }

}  // end namespace Chroma
