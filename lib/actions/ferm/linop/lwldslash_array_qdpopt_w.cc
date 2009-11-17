// $Id: lwldslash_array_qdpopt_w.cc,v 3.6 2009/11/14 20:01:46 eneil Exp $
/*! \file
 *  \brief Wilson Dslash linear operator array
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_array_qdpopt_w.h"


namespace Chroma 
{ 
  //! General Wilson-Dirac dslash of arrays
  /*! \ingroup linop
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */


  //! Creation routine
  void QDPWilsonDslashArrayOpt::create(Handle< FermState<T,P,Q> > state, int N5_)
  {
    multi1d<Real> cf(Nd);
    cf = 1.0;
    create(state, N5_, cf);
  }


  //! Creation routine with anisotropy
  void QDPWilsonDslashArrayOpt::create(Handle< FermState<T,P,Q> > state, int N5_,
				       const AnisoParam_t& anisoParam) 
  {
    START_CODE();

    create(state, N5_, makeFermCoeffs(anisoParam));

    END_CODE();
  }

  //! Creation routine
  void QDPWilsonDslashArrayOpt::create(Handle< FermState<T,P,Q> > state, int N5_,
				       const multi1d<Real>& coeffs_)
  {
    START_CODE();

    N5 = N5_;
    coeffs = coeffs_;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "QDPWilsonDslashArrayOpt: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Get links
    u = state->getLinks();

    // Rescale the u fields by the anisotropy
    for(int mu=0; mu < u.size(); ++mu)
    {
      u[mu] *= coeffs[mu];
    }

    END_CODE();
  }


  //! General Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi      Result				                (Write)
   *  \param psi      Pseudofermion field				(Read)
   *  \param isign    D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
   *  \param cb	      Checkerboard of OUTPUT vector			(Read) 
   */
  void 
  QDPWilsonDslashArrayOpt::apply (multi1d<LatticeFermion>& chi, 
			       const multi1d<LatticeFermion>& psi, 
			       enum PlusMinus isign, int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    for(int n=0; n < N5; ++n)
      apply(chi[n], psi[n], isign, cb);

    END_CODE();
  }



  //! General Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi	      Result				                (Write)
   *  \param psi	      Pseudofermion field				(Read)
   *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
   *  \param cb	      Checkerboard of OUTPUT vector			(Read) 
   */
  void 
  QDPWilsonDslashArrayOpt::apply (LatticeFermion& chi, const LatticeFermion& psi, 
			       enum PlusMinus isign, int cb) const
  {
    START_CODE();
#if (QDP_NC == 2) || (QDP_NC == 3)
    /*     F 
     *   a2  (x)  :=  U  (x) (1 - isign gamma  ) psi(x)
     *     mu          mu                    mu
     */
    /*     B           +
     *   a2  (x)  :=  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
     *     mu          mu                       mu
     */
    // Recontruct the bottom two spinor components from the top two
    /*                        F           B
     *   chi(x) :=  sum_mu  a2  (x)  +  a2  (x)
     *                        mu          mu
     */

    /* Why are these lines split? An array syntax would help, but the problem is deeper.
     * The expression templates require NO variable args (like int's) to a function
     * and all args must be known at compile time. Hence, the function names carry
     * (as functions usually do) the meaning (and implicit args) to a function.
     */
    int otherCB = (cb == 0 ? 1 : 0);

    switch (isign)
    {
    case PLUS: 
      {


	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;

	tmp[rb[otherCB]]  = spinProjectDir0Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 0);
	chi[rb[cb]] = spinReconstructDir0Minus(u[0]*tmp2);
	
	
	tmp[rb[otherCB]]  = spinProjectDir1Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Minus(u[1]*tmp2);

	tmp[rb[otherCB]]  = spinProjectDir2Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Minus(u[2]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir3Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Minus(u[3]*tmp2);
	
	
	tmp[rb[otherCB]]  = adj(u[0])*spinProjectDir0Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb[cb]] += spinReconstructDir0Plus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[1])*spinProjectDir1Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Plus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[2])*spinProjectDir2Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Plus(tmp2);

	tmp[rb[otherCB]]  = adj(u[3])*spinProjectDir3Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Plus(tmp2);	
	
      }


      break;

    case MINUS:
      {


	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;

	tmp[rb[otherCB]]  = spinProjectDir0Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 0);
	chi[rb[cb]] = spinReconstructDir0Plus(u[0]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir1Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Plus(u[1]*tmp2);

	tmp[rb[otherCB]] = spinProjectDir2Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Plus(u[2]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir3Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Plus(u[3]*tmp2);
	
	
	tmp[rb[otherCB]]  = adj(u[0])*spinProjectDir0Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb[cb]] += spinReconstructDir0Minus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[1])*spinProjectDir1Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Minus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[2])*spinProjectDir2Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Minus(tmp2);

	tmp[rb[otherCB]]  = adj(u[3])*spinProjectDir3Minus(psi);    
	tmp2 = shift(tmp, BACKWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Minus(tmp2);	
	
      }

      break;
    }

    getFermBC().modifyF(chi, QDP::rb[cb]);
#else
    QDPIO::cerr<<"lwldslash_array_qdpopt_w: not implemented for NC!=3\n";
    QDP_abort(13) ;
#endif

    END_CODE();
  }

} // End Namespace Chroma

