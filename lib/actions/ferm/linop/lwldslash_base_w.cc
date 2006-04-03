// $Id: lwldslash_base_w.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_base_w.h"
#include "io/aniso_io.h"


namespace Chroma 
{ 
  //! Take deriv of D
  /*!
   * \param chi     left vector                                 (Read)
   * \param psi     right vector                                (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void
  WilsonDslashBase::deriv(multi1d<LatticeColorMatrix>& ds_u,
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    deriv(ds_u, chi, psi, isign, 0);
    deriv(ds_tmp, chi, psi, isign, 1);
    ds_u += ds_tmp;

    END_CODE();
  }


  //! Take deriv of D
  /*! \return Computes   chi^dag * \dot(D} * psi  */
  void 
  WilsonDslashBase::deriv(multi1d<LatticeColorMatrix>& ds_u,
			  const LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign, int cb) const
  {
    START_CODE();

    ds_u.resize(Nd);

    AnisoParam_t anisoParam = getAnisoParam();
    multi1d<Real> anisoWeights(Nd);
    anisoWeights = 1;

    Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));

    if (anisoParam.anisoP)
    {
      // Set the weights
      for(int mu=0; mu < Nd; ++mu)
      {
	if (mu != anisoParam.t_dir)
	  anisoWeights[mu] *= ff;
      }
    }


    switch (isign)
    {
    case PLUS:
      for(int mu = 0; mu < Nd; ++mu)
      {
	// Undaggered:
        ds_u[mu][rb[cb]]   = anisoWeights[mu] * traceSpin(outerProduct(shift(psi - Gamma(1 << mu)*psi, FORWARD, mu),chi));
	ds_u[mu][rb[1-cb]] = zero;

	// The piece that comes from the U^daggered term. 
	// This piece is just -dagger() of the piece from applying
	// this function on the opposite checkerboard. It therefore
	// only contributes a factor of 2 to the traceless antihermitian
	// part of the result. Which should be swept into the taproj
	// normalisation. Right now until then, I explicitly multiply
	// the result by 0.5 below.

	// ds_u[mu][rb[1-cb]]  = traceSpin(outerProduct(psi + Gamma(1 << mu)*psi,shift(chi, FORWARD, mu)))*adj(u[mu]);
       	// ds_u[mu][rb[1-cb]] *= -Real(1);
	//
	// From factor of 2 that comes from the U^daggered piece.
	// This maybe should be absorbed into the taproj normalisation
	//
	// ds_u[mu] *= Real(0.5);
      }
      break;

    case MINUS:
      for(int mu = 0; mu < Nd; ++mu)
      {
	// Daggered:
	ds_u[mu][rb[cb]]   = anisoWeights[mu] * traceSpin(outerProduct(shift(psi + Gamma(1 << mu)*psi, FORWARD, mu),chi));
	
	ds_u[mu][rb[1-cb]] = zero;
	
	// The piece that comes from the U^daggered term. 
	// This piece is just -dagger() of the piece from applying
	// this function on the opposite checkerboard. It therefore
	// only contributes a factor of 2 to the traceless antihermitian
	// part of the result. Which should be swept into the taproj
	// normalisation. Right now until then, I explicitly multiply
	// the result by 0.5 below.
	//
	//        ds_u[mu][rb[1-cb]]  = traceSpin(outerProduct(psi - Gamma(1 << mu)*psi,shift(chi, FORWARD, mu)))*adj(u[mu]);
	//        ds_u[mu][rb[1-cb]] *= -Real(1);
	//	 
	// From factor of 2 that comes from the U^daggered piece.
	// This maybe should be absorbed into the taproj normalisation
	//
	// ds_u[mu] *= Real(0.5);
      }
      break;

    default:
      QDP_error_exit("unknown case");
    }
    
    getFermBC().zero(ds_u);

    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long 
  WilsonDslashBase::nFlops() const {return 1320;}

}; // End Namespace Chroma

