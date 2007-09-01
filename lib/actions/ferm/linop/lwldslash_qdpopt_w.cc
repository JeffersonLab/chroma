// $Id: lwldslash_qdpopt_w.cc,v 3.4 2007-09-01 23:44:10 uid3790 Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_qdpopt_w.h"


namespace Chroma 
{ 
  //! General Wilson-Dirac dslash
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

  //! Empty constructor
  QDPWilsonDslashOpt::QDPWilsonDslashOpt() {}
  
  //! Full constructor
  QDPWilsonDslashOpt::QDPWilsonDslashOpt(Handle< FermState<T,P,Q> > state)
  {
    create(state);
  }
  
  //! Full constructor with anisotropy
  QDPWilsonDslashOpt::QDPWilsonDslashOpt(Handle< FermState<T,P,Q> > state,
				   const AnisoParam_t& aniso_) 
  {
    create(state, aniso_);
  }

  //! Creation routine
  void QDPWilsonDslashOpt::create(Handle< FermState<T,P,Q> > state)
  {
    AnisoParam_t foo;
    create(state, foo);
  }

  //! Creation routine with anisotropy
  void QDPWilsonDslashOpt::create(Handle< FermState<T,P,Q> > state,
			       const AnisoParam_t& aniso_) 
  {
    anisoParam = aniso_;
    fbc = state->getFermBC();
    u   = state->getLinks();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "WilsonDslash: error: fbc is null" << endl;
      QDP_abort(1);
    }

    if (anisoParam.anisoP)
    {
      Real ff = anisoParam.nu / anisoParam.xi_0;
  
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != anisoParam.t_dir)
	  u[mu] *= ff;
      }
    }

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
  QDPWilsonDslashOpt::apply (LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign, int cb) const
  {
    START_CODE();

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

    END_CODE();
  }

} // End Namespace Chroma

