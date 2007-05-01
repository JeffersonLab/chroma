// $Id: invmr.cc,v 3.4 2007-05-01 15:27:27 bjoo Exp $
/*! \file
 *  \brief Minimal-Residual (MR) for a generic fermion Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invmr.h"

using namespace QDP::Hints;

namespace Chroma 
{

  //! Minimal-residual (MR) algorithm for a generic Linear Operator 
  /*! \ingroup invert
   * This subroutine uses the Minimal Residual (MR) algorithm to determine
   * the solution of the set of linear equations. Here we allow M to be nonhermitian.
   *
   *   	    Chi  =  M . Psi 
   *
   * Algorithm:
   *
   *  Psi[0]                                      Argument
   *  r[0]    :=  Chi  -  M . Psi[0] ;            Initial residual
   *  IF |r[0]| <= RsdCG |Chi| THEN RETURN;       Converged?
   *  FOR k FROM 1 TO MaxCG DO                    MR iterations
   *      a[k-1]  := <M.r[k-1],r[k-1]> / <M.r[k-1],M.r[k-1]> ;
   *      ap[k-1] := MRovpar * a[k] ;             Overrelaxtion step
   *      Psi[k]  += ap[k-1] r[k-1] ;   	        New solution vector
   *      r[k]    -= ap[k-1] A . r[k-1] ;         New residual
   *      IF |r[k]| <= RsdCG |Chi| THEN RETURN;   Converged?

   * Arguments:

   *  \param M       Linear Operator             (Read)
   *  \param chi     Source                      (Read)
   *  \param psi     Solution                    (Modify)
   *  \param RsdCG   MR residual accuracy        (Read)
   *  \param MRovpar Overrelaxation parameter    (Read)
   *  \param MaxMR   Maximum MR iterations       (Read)

   * Local Variables:

   *  r   	Residual vector
   *  cp  	| r[k] |**2
   *  c   	| r[k-1] |**2
   *  k   	MR iteration counter
   *  a   	a[k]
   *  d   	< M.r[k], M.r[k] >
   *  R_Aux     Temporary for  M.Psi
   *  Mr        Temporary for  M.r

   * Global Variables:

   *  MaxMR       Maximum number of MR iterations allowed
   *  RsdCG       Maximum acceptable MR residual (relative to source)
   *
   * Subroutines:
   *
   *  M           Apply matrix to vector
   *
   * @{
   */

  template<typename T, typename C>
  SystemSolverResults_t 
  InvMR_a(const C& M,
	  const T& chi,
	  T& psi,
	  const Real& MRovpar,
	  const Real& RsdMR, 
	  int MaxMR, 
	  enum PlusMinus isign)
  {
    START_CODE();

    const Subset& s = M.subset();

    SystemSolverResults_t  res;
    moveToFastMemoryHint(psi,true);
    T Mr;                moveToFastMemoryHint(Mr);
    T chi_internal;      moveToFastMemoryHint(chi_internal);

    chi_internal[s] = chi;

    Complex a;
    DComplex c;
    Double d;
    int k;

    QDPIO::cout << "InvMR: starting" << endl;
    FlopCounter flopcount;
    flopcount.reset();
    StopWatch swatch;
    swatch.reset();
    swatch.start();

    Real rsd_sq = (RsdMR * RsdMR) * Real(norm2(chi_internal,s));
    flopcount.addSiteFlops(4*Nc*Ns,s);
        
    /*  r[0]  :=  Chi - M . Psi[0] */
    /*  r  :=  M . Psi  */
    M(Mr, psi, isign);
    flopcount.addFlops(M.nFlops());

    T r;                moveToFastMemoryHint(r);
    r[s] = chi_internal - Mr;
    flopcount.addSiteFlops(2*Nc*Ns,s);
  
    /*  Cp = |r[0]|^2 */
    Double cp = norm2(r, s);                 /* 2 Nc Ns  flops */
    flopcount.addSiteFlops(4*Nc*Ns, s);

//  QDPIO::cout << "InvMR: k = 0  cp = " << cp << "  rsd_sq = " << rsd_sq << endl;

    /*  IF |r[0]| <= RsdMR |Chi| THEN RETURN; */
    if ( toBool(cp  <=  rsd_sq) )
    {
      res.n_count = 0;
      res.resid   = sqrt(cp);
      swatch.stop();
      flopcount.report("invMR", swatch.getTimeInSeconds());
      revertFromFastMemoryHint(psi,true);
      END_CODE();
      return res;
    }

    /*  FOR k FROM 1 TO MaxMR DO */
    k = 0;
    while( (k < MaxMR) && (toBool(cp > rsd_sq)) )
    {
      ++k;

      /*  a[k-1] := < M.r[k-1], r[k-1] >/ < M.r[k-1], M.r[k-1] > ; */
      /*  Mr = M * r  */
      M(Mr, r, isign);  flopcount.addFlops(M.nFlops());

      /*  c = < M.r, r > */
      c = innerProduct(Mr, r, s);  flopcount.addSiteFlops(4*Nc*Ns,s);
    
      /*  d = | M.r | ** 2  */
      d = norm2(Mr, s);  flopcount.addSiteFlops(4*Nc*Ns,s);

      /*  a = c / d */
      a = c / d;
    
      /*  a[k-1] *= MRovpar ; */
      a = a * MRovpar;

      /*  Psi[k] += a[k-1] r[k-1] ; */
      psi[s] += r * a;    flopcount.addSiteFlops(4*Nc*Ns,s);

      /*  r[k] -= a[k-1] M . r[k-1] ; */
      r[s] -= Mr * a;    flopcount.addSiteFlops(4*Nc*Ns,s);

      /*  cp  =  | r[k] |**2 */
      cp = norm2(r, s);    flopcount.addSiteFlops(4*Nc*Ns,s);

//    QDPIO::cout << "InvMR: k = " << k << "  cp = " << cp << endl;
    }
    res.n_count = k;
    res.resid   = sqrt(cp);
    swatch.stop();
    QDPIO::cout << "InvMR: k = " << k << "  cp = " << cp << endl;
    flopcount.report("invmr", swatch.getTimeInSeconds());
    revertFromFastMemoryHint(psi,true);

    // Compute the actual residual
    {
      M(Mr, psi, isign);
      Double actual_res = norm2(chi_internal - Mr,s);
      res.resid = sqrt(actual_res);
    }

    if ( res.n_count == MaxMR )
      QDPIO::cerr << "Nonconvergence Warning" << endl;
    
    END_CODE();
    return res;
  }


  // Fix here for now
  template<>
  SystemSolverResults_t 
  InvMR(const LinearOperator<LatticeFermion>& M,
	const LatticeFermion& chi,
	LatticeFermion& psi,
	const Real& MRovpar,
	const Real& RsdMR, 
	int MaxMR,
	enum PlusMinus isign)
  {
    return InvMR_a(M, chi, psi, MRovpar, RsdMR, MaxMR, isign);
  }


  // Fix here for now
  template<>
  SystemSolverResults_t 
  InvMR(const LinearOperator<LatticeStaggeredFermion>& M,
	const LatticeStaggeredFermion& chi,
	LatticeStaggeredFermion& psi,
	const Real& MRovpar,
	const Real& RsdMR, 
	int MaxMR,
	enum PlusMinus isign)
  {
    return InvMR_a(M, chi, psi, MRovpar, RsdMR, MaxMR, isign);
  }

  /*! @} */  // end of group invert

}  // end namespace Chroma
