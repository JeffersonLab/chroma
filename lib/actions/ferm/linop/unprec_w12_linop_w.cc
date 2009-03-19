// $Id: unprec_w12_linop_w.cc,v 3.2 2009-03-19 17:08:00 mcneile Exp $
/*! \file
 *  \brief Unpreconditioned W12 action
 */

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3


#include "chromabase.h"
#include "actions/ferm/linop/unprec_w12_linop_w.h"

namespace Chroma 
{ 
  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  parameters   	       (Read)
   */
  void UnprecW12LinOp::create(Handle< FermState<T,P,Q> > fs,
			      const CloverFermActParams& param_)
  {
    u = fs->getLinks();
    param = param_;

    A.create(fs, param); // Unaltered fields for clover

    if (Nd != 4)
    {
      QDPIO::cerr << "UnprecW12 only supports Nd=4" << endl;
      QDP_abort(1);
    }

    Real ff;
  
    if (param.anisoParam.anisoP)
    {
      aniso_fact = param.anisoParam.nu / param.anisoParam.xi_0;
      j_decay = param.anisoParam.t_dir;
    }
    else
    {
      aniso_fact = 1.0;
      j_decay = Nd;
    }

    fact1 = 1 + (Nd-1)*aniso_fact + param.Mass;
    fact2 = -2 /(3*param.u0);
    fact4 = 1 / (12*param.u0*param.u0); 
    fact3 = 2 * fact4;

    if (j_decay != Nd-1)
    {
      QDPIO::cerr << "W12: implementation restriction: j_decay must be Nd-1" << endl;
      QDP_abort(1);
    }

    // Anisotropy not yet folded into gauge fields
  }


  //! GAMW
  /*!
   * Description:
   *
   * This routine applies the operator W' to Psi, putting the result in Chi.
   *
   *   chi(x,mu)  :=  + (1 - isign gamma  ) U  (x) psi(x+mu)
   *                                    mu   mu
   *                                         +
   *                  + (1 + isign gamma  ) U  (x-mu) psi(x-mu)
   *                                    mu   mu
   */
  void UnprecW12LinOp::gamW(multi1d<LatticeFermion>& chi, 
			    const LatticeFermion& psi,
			    int j_decay,
			    enum PlusMinus isign) const
  {
    START_CODE();
  
    chi.resize(Nd);

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
    switch (isign)
    {
    case PLUS:
      chi[0] = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi), FORWARD, 0))
  	     + spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi), BACKWARD, 0));
      chi[1] = spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi), FORWARD, 1))
	     + spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi), BACKWARD, 1));
      chi[2] = spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi), FORWARD, 2))
	     + spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi), BACKWARD, 2));
      chi[3] = spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi), FORWARD, 3))
	     + spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi), BACKWARD, 3));
      break;

    case MINUS:
      chi[0] = spinReconstructDir0Plus(u[0] * shift(spinProjectDir0Plus(psi), FORWARD, 0))
	     + spinReconstructDir0Minus(shift(adj(u[0]) * spinProjectDir0Minus(psi), BACKWARD, 0));
      chi[1] = spinReconstructDir1Plus(u[1] * shift(spinProjectDir1Plus(psi), FORWARD, 1))
	     + spinReconstructDir1Minus(shift(adj(u[1]) * spinProjectDir1Minus(psi), BACKWARD, 1));
      chi[2] = spinReconstructDir2Plus(u[2] * shift(spinProjectDir2Plus(psi), FORWARD, 2))
	     + spinReconstructDir2Minus(shift(adj(u[2]) * spinProjectDir2Minus(psi), BACKWARD, 2));
      chi[3] = spinReconstructDir3Plus(u[3] * shift(spinProjectDir3Plus(psi), FORWARD, 3))
	     + spinReconstructDir3Minus(shift(adj(u[3]) * spinProjectDir3Minus(psi), BACKWARD, 3));
      break;
    }
      
    END_CODE();
  }


  //! GAMWMU
  /*! This routine applies the operator W' to Psi, putting the result in Chi.
   *
   *   chi(x,mu)  :=  + (1 - isign gamma  ) U  (x) psi  (x+mu)
   *                                    mu   mu       mu
   *                                         +
   *                  + (1 + isign gamma  ) U  (x-mu) psi  (x-mu)
   *                                    mu   mu          mu
   */
  void UnprecW12LinOp::gamWmu(multi1d<LatticeFermion>& chi, 
			      const multi1d<LatticeFermion>& psi,
			      int j_decay,
			      enum PlusMinus isign) const
  {
    START_CODE();
  
    chi.resize(Nd);

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
    switch (isign)
    {
    case PLUS:
      chi[0] = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi[0]), FORWARD, 0))
  	     + spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi[0]), BACKWARD, 0));
      chi[1] = spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi[1]), FORWARD, 1))
	     + spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi[1]), BACKWARD, 1));
      chi[2] = spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi[2]), FORWARD, 2))
	     + spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi[2]), BACKWARD, 2));
      chi[3] = spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi[3]), FORWARD, 3))
	     + spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi[3]), BACKWARD, 3));
      break;

    case MINUS:
      chi[0] = spinReconstructDir0Plus(u[0] * shift(spinProjectDir0Plus(psi[0]), FORWARD, 0))
	     + spinReconstructDir0Minus(shift(adj(u[0]) * spinProjectDir0Minus(psi[0]), BACKWARD, 0));
      chi[1] = spinReconstructDir1Plus(u[1] * shift(spinProjectDir1Plus(psi[1]), FORWARD, 1))
	     + spinReconstructDir1Minus(shift(adj(u[1]) * spinProjectDir1Minus(psi[1]), BACKWARD, 1));
      chi[2] = spinReconstructDir2Plus(u[2] * shift(spinProjectDir2Plus(psi[2]), FORWARD, 2))
	     + spinReconstructDir2Minus(shift(adj(u[2]) * spinProjectDir2Minus(psi[2]), BACKWARD, 2));
      chi[3] = spinReconstructDir3Plus(u[3] * shift(spinProjectDir3Plus(psi[3]), FORWARD, 3))
	     + spinReconstructDir3Minus(shift(adj(u[3]) * spinProjectDir3Minus(psi[3]), BACKWARD, 3));
      break;
    }
      
    END_CODE();
  }


  // Override inherited one with a few more funkies
  /*!
   * Chi = (m0 - (2/3*((1/2)*(1/4))*sigma.F + W'  + (1/6)*W^2_mu) * Psi
   *
   * For sake of generality, efficiency, and future tinkering, this routine
   * actually does
   *
   *  Chi =  (clov + Kappa_ds*W'  + Kappa_w2*W^2_mu) * Psi    
   *
   * where  
   *
   *   clov = 1 + const*sigma.F
   *
   * and const is some constant depending on kappa described in  conslinop
   */
  void UnprecW12LinOp::operator()(LatticeFermion & chi, 
				  const LatticeFermion& psi,
				  enum PlusMinus isign) const
  {
    START_CODE();

    multi1d<LatticeFermion> tmp1(Nd);
    multi1d<LatticeFermion> tmp2(Nd);

    // Apply clover term
    A(chi, psi, isign);

    // Add the diagonal term from W^2
    chi += fact1 * psi;

    /*  Tmp1_mu  =  D'_mu * Psi */
    /*  Chi     -=  sum_mu Kappa_ds * Tmp1_mu */
    gamW(tmp1, psi, Nd, isign);
    for(int mu = 0; mu < Nd; ++mu)
    {
      if (mu != j_decay)
	chi -= Real(aniso_fact * fact2) * tmp1[mu];
      else
	chi -= fact2 * tmp1[mu];
    }

    /*  Tmp2_mu  =  D'_mu * Tmp1_mu */
    /*  Chi     +=  sum_mu  Kappa_w2 * Tmp1_mu */
    gamWmu(tmp2, tmp1, j_decay, isign);
    for(int mu = 0; mu < Nd; ++mu)
      if (mu != j_decay)
	chi += fact4 * tmp2[mu];

    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long UnprecW12LinOp::nFlops() const
  { 
    unsigned long site_flops = 0;
    return site_flops*(Layout::sitesOnNode());
  }


  //! Derivative of unpreconditioned W12 dM/dU
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  void 
  UnprecW12LinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			const LatticeFermion& chi, const LatticeFermion& psi, 
			enum PlusMinus isign) const
  {
    QDP_error_exit("W12 deriv not correct yet");


    START_CODE();

    // This does both parities
    A.deriv(ds_u, chi, psi, isign);

    // Factor in front of the dslash
    for(int mu = 0; mu < Nd; ++mu)
      ds_u[mu] *= fact2;

    // NOTE: missing derivative of 2 link piece

    END_CODE();
  }

}; // End Namespace Chroma


#endif
#endif
#endif
