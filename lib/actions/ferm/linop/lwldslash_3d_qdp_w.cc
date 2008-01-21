// $Id: lwldslash_3d_qdp_w.cc,v 3.2 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_3d_qdp_w.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

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
  QDPWilsonDslash3D::QDPWilsonDslash3D() {}
  
  //! Full constructor
  QDPWilsonDslash3D::QDPWilsonDslash3D(Handle< FermState<T,P,Q> > state)
  {
    create(state);
  }

  
  //! Full constructor with anisotropy
  QDPWilsonDslash3D::QDPWilsonDslash3D(Handle< FermState<T,P,Q> > state,
				       const AnisoParam_t& aniso_) 
  {
    create(state, aniso_);
  }

  //! Creation routine
  void QDPWilsonDslash3D::create(Handle< FermState<T,P,Q> > state)
  {
    multi1d<Real> cf(Nd);
    cf = 1.0;
    create(state, cf);
  }

  //! Creation routine with anisotropy
  void QDPWilsonDslash3D::create(Handle< FermState<T,P,Q> > state,
				 const AnisoParam_t& anisoParam) 
  {
    START_CODE();

    create(state, makeFermCoeffs(anisoParam));

    END_CODE();
  }

  //! Full constructor with general coefficients
  void QDPWilsonDslash3D::create(Handle< FermState<T,P,Q> > state,
				 const multi1d<Real>& coeffs_)
  {
    coeffs = coeffs_;
    fbc = state->getFermBC();
    u   = state->getLinks();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "WilsonDslash3D: error: fbc is null" << endl;
      QDP_abort(1);
    }
  
    // Rescale the u fields by the anisotropy
    for(int mu=0; mu < u.size(); ++mu)
    {
      u[mu] *= coeffs[mu];
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
  QDPWilsonDslash3D::apply (LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign, int cb) const
  {
    START_CODE();
    

    int otherCB = (cb == 0 ? 1 : 0);

    switch (isign)
    {
    case PLUS: 
      {


	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;

	tmp[rb3[otherCB]]  = spinProjectDir0Minus(psi);
	tmp2[rb3[cb]] = shift(tmp, FORWARD, 0);
	chi[rb3[cb]] = spinReconstructDir0Minus(u[0]*tmp2);
	
	
	tmp[rb3[otherCB]]  = spinProjectDir1Minus(psi);
	tmp2[rb3[cb]] = shift(tmp, FORWARD, 1);
	chi[rb3[cb]] += spinReconstructDir1Minus(u[1]*tmp2);

	tmp[rb3[otherCB]]  = spinProjectDir2Minus(psi);
	tmp2[rb3[cb]] = shift(tmp, FORWARD, 2);
	chi[rb3[cb]] += spinReconstructDir2Minus(u[2]*tmp2);
	
	
	
	tmp[rb3[otherCB]]  = adj(u[0])*spinProjectDir0Plus(psi);
	tmp2[rb3[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb3[cb]] += spinReconstructDir0Plus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[1])*spinProjectDir1Plus(psi);
	tmp2[rb3[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb3[cb]] += spinReconstructDir1Plus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[2])*spinProjectDir2Plus(psi);
	tmp2[rb3[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb3[cb]] += spinReconstructDir2Plus(tmp2);

	
      }


      break;

    case MINUS:
      {


	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;

	tmp[rb3[otherCB]]  = spinProjectDir0Plus(psi);
	tmp2[rb3[cb]] = shift(tmp, FORWARD, 0);
	chi[rb3[cb]] = spinReconstructDir0Plus(u[0]*tmp2);
	
	tmp[rb3[otherCB]]  = spinProjectDir1Plus(psi);
	tmp2[rb3[cb]] = shift(tmp, FORWARD, 1);
	chi[rb3[cb]] += spinReconstructDir1Plus(u[1]*tmp2);

	tmp[rb3[otherCB]] = spinProjectDir2Plus(psi);
	tmp2[rb3[cb]] = shift(tmp, FORWARD, 2);
	chi[rb3[cb]] += spinReconstructDir2Plus(u[2]*tmp2);
	
	
	tmp[rb3[otherCB]]  = adj(u[0])*spinProjectDir0Minus(psi);
	tmp2[rb3[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb3[cb]] += spinReconstructDir0Minus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[1])*spinProjectDir1Minus(psi);
	tmp2[rb3[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb3[cb]] += spinReconstructDir1Minus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[2])*spinProjectDir2Minus(psi);
	tmp2[rb3[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb3[cb]] += spinReconstructDir2Minus(tmp2);

	
      }

      break;
    }

    END_CODE();
  }

} // End Namespace Chroma

#endif
#endif
#endif
