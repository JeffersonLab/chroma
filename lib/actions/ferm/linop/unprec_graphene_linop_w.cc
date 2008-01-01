// $Id: unprec_graphene_linop_w.cc,v 1.4 2008-01-01 22:45:31 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Graphene fermion linear operator.
 *
 * This formulation follows Borici's variant of Creutz's graphene
 * fermion construction. Borici's variant is described in
 * arXiv:0712.4401 and Cruetz's original construction is described
 * in arXiv:0712.1201
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_graphene_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   */
  void UnprecGrapheneLinOp::create(Handle< FermState<T,P,Q> > fs,
				   const Real& Mass_)
  {
    AnisoParam_t anisoParam;
    create(fs, Mass_, anisoParam);
  }


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   * \param aniso   anisotropy struct   	       (Read)
   */
  void UnprecGrapheneLinOp::create(Handle< FermState<T,P,Q> > fs,
				   const Real& Mass_,
				   const AnisoParam_t& aniso_)
  {
    START_CODE();

    Mass = Mass_;
    anisoParam = aniso_;
    fbc = fs->getFermBC();
    u   = fs->getLinks();

    // Sanity checks
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "GrapheneLinOp: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Insist on this for graphene
    if (QDP::Nd != 4 || QDP::Ns != 4)
    {
      QDPIO::cerr << "GrapheneLinOp: requires Nd=4 and Ns=4" << endl;
      QDP_abort(1);
     
    }

    // For the moment, I'm hesitant to turn on anisotropy. I think the
    // construction of graphene goes through in this case. The point
    // is that the *fermion* anisotropy (really gamma_f=xi_0/nu in my
    // new language) is what is tuned so the spatial and temporal lattice
    // spacing are equal once the anisotropy factor is restored. However,
    // this all could stand some more thought, so disable it here.
    if (anisoParam.anisoP) 
    {
      QDPIO::cerr << "UnprecGraphene: not fully supporting anisotropy at this moment" << endl;
      QDPIO::cerr << "UnprecGraphene: only this one check is stopping it from functioning; otherwise, the code is in place." << endl;
      QDP_abort(1);
    }

    Real ff;
    if (anisoParam.anisoP) 
      ff = anisoParam.nu / anisoParam.xi_0;
    else
      ff = Real(1);
  
    if (anisoParam.anisoP)
    {
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != anisoParam.t_dir)
	  u[mu] *= ff;
      }
    }

    // This is Borici's "alpha" matrix. Here, Nd must be 4.
    alpha.resize(Nd,Nd);

    for(int mu=0; mu < Nd; ++mu)
    {
      for(int nu=0; nu < Nd; ++nu)
      {
	alpha(mu,nu) = (mu == nu) ? 1 : -1;
      }
    }
	
    END_CODE();
  }


  // Form  gamma_mu * psi
  void UnprecGrapheneLinOp::gammaMults(multi1d<LatticeFermion>& gams, 
				       const LatticeFermion& psi) const
  {
    gams.resize(Nd);   moveToFastMemoryHint(gams);

    // Build gamma matrix multiplied pieces
    gams[0] = GammaConst<Ns,1>()*psi;
    gams[1] = GammaConst<Ns,2>()*psi;
    gams[2] = GammaConst<Ns,4>()*psi;
    gams[3] = GammaConst<Ns,8>()*psi;
  }


  // Hop terms
  void UnprecGrapheneLinOp::iGamMu(LatticeFermion& iGam,
				   const multi1d<LatticeFermion>& gams,
				   int mu) const
  {
    LatticeFermion tmp2;   moveToFastMemoryHint(tmp2);

    // The Gamma piece. This will be shift later. Unroll the loop.
    // Flops could be saved here since some pieces are re-added again
    // because of the similar sign structure of alpha.
    tmp2 = Real(alpha(mu,0))*gams[0];
    for(int nu=1; nu < Nd; ++nu)
      tmp2 += Real(alpha(mu,nu))*gams[nu];

    // i*Gamma_mu * psi
    iGam = timesI(tmp2);
  }

  
  //! Apply unpreconditioned Graphene fermion linear operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field     	       (Read)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecGrapheneLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
					enum PlusMinus isign) const
  {
    START_CODE();

    //
    //  Chi   =  (Nd+Mass)*Psi  -  (1/2) * D' Psi
    //
    multi1d<LatticeFermion> gams;
    LatticeFermion tmp2;  moveToFastMemoryHint(tmp2);
    LatticeFermion tmp3;  moveToFastMemoryHint(tmp3);
    LatticeFermion iGam;  moveToFastMemoryHint(iGam);
    Real half = 0.5;

    // Build gamma matrix multiplied pieces
    gammaMults(gams, psi);

    // Mass term piece. Unroll.
    tmp2 = gams[0];
    for(int mu=1; mu < Nd; ++mu)
      tmp2 += gams[mu];

    chi = timesI(tmp2);

    // Hop pieces
    for(int mu=0; mu < Nd; ++mu)
    {
      // Unshifted hop terms
      iGamMu(iGam, gams, mu);

      // Forward piece. 
      tmp2 = iGam + gams[mu];
      tmp3 = u[mu] * shift(tmp2, FORWARD, mu);
      chi += half * tmp3;

      // Backward piece.
      tmp2 = iGam - gams[mu];
      tmp3 = adj(u[mu]) * tmp2;
      chi += half * shift(tmp3, BACKWARD, mu);
    }
    
    if (isign == MINUS)
      chi = -chi;

    chi += Mass*psi;

    getFermBC().modifyF(chi);
  
    END_CODE();
  }



  //! Derivative of unpreconditioned Graphene dM/dU
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  void 
  UnprecGrapheneLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			     const LatticeFermion& chi, const LatticeFermion& psi, 
			     enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    // Fold in the 1/2 from the front of the hop terms here.
    multi1d<Real> anisoWeights(Nd);
    anisoWeights = 0.5;

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

    // Build gamma matrix multiplied pieces
    multi1d<LatticeFermion> gams;
    LatticeFermion tmp2;  moveToFastMemoryHint(tmp2);
    LatticeFermion iGam;  moveToFastMemoryHint(iGam);
 
    gammaMults(gams, psi);

    for(int mu = 0; mu < Nd; ++mu) 
    {
      // Unshifted hop terms
      iGamMu(iGam, gams, mu);

      switch (isign) 
      {
      case PLUS:
      {
	// Forward piece. 
	tmp2 = iGam + gams[mu];
      }
      break;

      case MINUS:
      {
	// Backward piece.
	// NOTE: This is confusing. There is no derivative of the U^dag, just the U.
	// However, we take D^dag here so this is gamma_5*D*gamma_5. For the non-mass
	// term, this is just -D. So, the derivative still picks up just for the forward
	// term, but with a minus sign.
	tmp2 = iGam + gams[mu];
	tmp2 = -tmp2;
      }
      break;

      default:
	QDP_error_exit("unknown case");
      }

      // This is the forward piece
      LatticeFermion tmp3 = shift(tmp2, FORWARD, mu);
      
      // This step supposedly optimised in QDP++
      LatticeColorMatrix temp_mat = traceSpin(outerProduct(tmp3,chi));
    
      // Just do the bit we need.
      ds_u[mu] = anisoWeights[mu] * temp_mat;
    }

    getFermBC().zero(ds_u);

    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long UnprecGrapheneLinOp::nFlops() const
  {
    // This is not quite right, but I don't think it's miles off.
    unsigned long site_flops = 4*Nc*Ns + 2*(10*Nc*Ns+8*264);
    return site_flops*Layout::sitesOnNode();
  }

} // End Namespace Chroma
