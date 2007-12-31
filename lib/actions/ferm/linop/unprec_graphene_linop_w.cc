// $Id: unprec_graphene_linop_w.cc,v 1.1 2007-12-31 23:24:26 edwards Exp $
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
    multi1d<LatticeFermion> tmp1(Nd);   moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2;                moveToFastMemoryHint(tmp2);
    LatticeFermion tmp3;                moveToFastMemoryHint(tmp3);
    LatticeFermion tmp4;                moveToFastMemoryHint(tmp4);
    Real half = 0.5;

    chi = Mass*psi;
 
    // Build gamma matrix multiplied pieces
    tmp1[0] = GammaConst<Ns,1>()*psi;
    tmp1[1] = GammaConst<Ns,2>()*psi;
    tmp1[2] = GammaConst<Ns,4>()*psi;
    tmp1[3] = GammaConst<Ns,8>()*psi;

    // Mass term piece. Unroll.
    tmp2 = tmp1[0];
    for(int mu=1; mu < Nd; ++mu)
      tmp2 += tmp1[mu];

    chi += timesI(tmp2);

    // Hop pieces
    for(int mu=0; mu < Nd; ++mu)
    {
      // The Gamma piece. This will be shift later. Unroll the loop.
      // Flops could be saved here since some pieces are re-added again
      // because of the similar sign structure of alpha.
      tmp2 = Real(alpha(mu,0))*tmp1[0];
      for(int nu=1; nu < Nd; ++nu)
	tmp2 += Real(alpha(mu,nu))*tmp1[nu];

      // NOTE: multiply by half here to keep multiply
      // unit busy, otherwise it would sit idle till another pass for
      // the multiply to be used.
      // Forward piece. 
      tmp3 = timesI(tmp2) + tmp1[mu];
      tmp4 = u[mu] * shift(tmp3, FORWARD, mu);
      chi += half * tmp4;

      // Backward piece.
      tmp3 = timesI(tmp2) - tmp1[mu];
      tmp4 = adj(u[mu]) * tmp3;
      chi += half * shift(tmp4, BACKWARD, mu);
    }
    
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

    QDPIO::cerr << "UnprecGraphene: deriv not implemented" << endl;
    QDP_abort(1);

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
