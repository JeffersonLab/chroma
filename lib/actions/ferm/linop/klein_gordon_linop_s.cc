// $Id: klein_gordon_linop_s.cc,v 1.3 2006-12-09 22:22:07 edwards Exp $
/*! \file
 *  \brief Klein-Gordon boson action masquerading action as a staggered action
 */

#include "chromabase.h"
#include "actions/ferm/linop/klein_gordon_linop_s.h"

namespace Chroma 
{ 
  // Full constructor
  KleinGordonLinOp::KleinGordonLinOp(Handle< FermState<T,P,Q> > state_,
				     const Real& Mass_,
				     const AnisoParam_t& aniso)
  {
    create(state_, Mass_, aniso);
  }


  // Creation routine with anisotropy
  void KleinGordonLinOp::create(Handle< FermState<T,P,Q> > state,
				const Real& Mass_,
				const AnisoParam_t& aniso_) 
  {
    START_CODE();

    // Save a copy of the aniso params original fields and with aniso folded in
    Mass = Mass_;
    anisoParam = aniso_;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "KleinGordonLinOp: error - fbc is null" << endl;
      QDP_abort(1);
    }

    // Fold in anisotropy
    u = state->getLinks();
    Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
    fact = 1 + (Nd-1)*ff + Mass;
  
    if (anisoParam.anisoP)
    {
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != anisoParam.t_dir)
	  u[mu] *= ff;
      }
    }

    END_CODE();
  }


  // Apply unpreconditioned Klein-Gordon operator
  void KleinGordonLinOp::operator()(LatticeStaggeredFermion& chi, 
				    const LatticeStaggeredFermion& psi, 
				    enum PlusMinus isign) const
  {
    START_CODE();

    LatticeStaggeredFermion tmp;   moveToFastMemoryHint(tmp);
    tmp = zero;

    for(int mu = 0; mu < Nd; ++mu )
    {
      tmp += u[mu]*shift(psi, FORWARD, mu);
      tmp += shift(adj(u[mu])*psi, BACKWARD, mu);
    }

    Real mhalf = -0.5;
    chi = fact*psi + mhalf*tmp;

    getFermBC().modifyF(chi);
  
    END_CODE();
  }


  // Derivative of unpreconditioned Klein-Gordon operator
  void 
  KleinGordonLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			  const LatticeStaggeredFermion& chi, 
			  const LatticeStaggeredFermion& psi, 
			  enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    multi1d<Real> anisoWeights(Nd);
    anisoWeights = 1;
    Real mhalf = -0.5;

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


    // The normalizations here are confusing (at least to RGE). See the comments
    // in  lwldslash_base_w.cc
    switch (isign)
    {
    case PLUS:
      for(int mu = 0; mu < Nd; ++mu)
      {
	// Undaggered:
        ds_u[mu] = mhalf * anisoWeights[mu] * traceSpin(outerProduct(shift(psi, FORWARD, mu),chi));
      }
      break;

    case MINUS:
      for(int mu = 0; mu < Nd; ++mu)
      {
	// Daggered:
	ds_u[mu] = mhalf * anisoWeights[mu] * traceSpin(outerProduct(shift(psi, FORWARD, mu),chi));
      }
      break;

    default:
      QDP_error_exit("unknown case");
    }
    
    getFermBC().zero(ds_u);

    END_CODE();
  }


  // Return flops performed by the operator()
  unsigned long KleinGordonLinOp::nFlops() const
  {
    // fact*psi + U*psi = (fact*psi = 6M) + (8 dirs)*(3 rows)*(3 cols)*(2M + 2A)
    unsigned long site_flops = 294;
    return site_flops*(Layout::sitesOnNode());
  }

} // End Namespace Chroma

