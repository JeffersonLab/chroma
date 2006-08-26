// $Id: unprec_hamberwu_linop_w.cc,v 3.1 2006-08-26 05:50:06 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Hamber-Wu linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_hamberwu_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 

  //! Creation routine
  /*!
   *
   * \param u_ 	    gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   */
  void UnprecHamberWuLinOp::create(Handle< FermState<T,P,Q> > fs,
				   const Real& Mass_, const Real& u0_)
  {
    Mass = Mass_;
    u0   = u0_;

    D.create(fs);

    const multi1d<LatticeColorMatrix>& u = fs->getLinks();

    //    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;

    /* Make "dble links", UU */
    u_dble.resize(u.size());
    for(int mu=0; mu < u.size(); ++mu)
      u_dble[mu] = u[mu] * shift(u[mu], FORWARD, mu);

    fact1 = Nd + Mass;
    fact2 = -2 /(3*u0);
    fact4 = 1 / (12*u0*u0); 
    fact3 = 2 * fact4;
  }



  // Override inherited one with a few more funkies
  void UnprecHamberWuLinOp::operator()(LatticeFermion & chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    LatticeFermion tmp1;    moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2;    moveToFastMemoryHint(tmp2);
    LatticeFermion tmp_s;   moveToFastMemoryHint(tmp_s);
    LatticeFermion tmp_v;   moveToFastMemoryHint(tmp_v);

    D(tmp1, psi, isign);

    chi = fact1*psi + fact2*tmp1;

    tmp_s = zero;
    tmp_v = zero;

    for (int mu = 0; mu < Nd; mu++)
    {
      /* Next-nearest neighbors */
      tmp2 = shift(psi, FORWARD, mu);
      tmp1 = u_dble[mu] * shift(tmp2, FORWARD, mu);
      tmp_s += tmp1;
      tmp2 = Gamma(1 << mu) * tmp1;
      tmp_v -= tmp2;

      tmp2 = shift(adj(u_dble[mu]) * psi, BACKWARD, mu);
      tmp1 = shift(tmp2, BACKWARD, mu);
      tmp_s += tmp1;
      tmp2 = Gamma(1 << mu) * tmp1;
      tmp_v += tmp2;
    }

    chi += fact3 * tmp_s;
    if (isign == PLUS)
      chi += fact4 * tmp_v;
    else
      chi -= fact4 * tmp_v;

    getFermBC().modifyF(chi);

    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long UnprecHamberWuLinOp::nFlops() const
  { 
    unsigned long site_flops = D.nFlops()+2*(2*Nc)*(2*Nc)*Ns*4*Nd;
    return site_flops*(Layout::sitesOnNode());
  }


  //! Derivative of unpreconditioned Hamber-Wu dM/dU
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  void 
  UnprecHamberWuLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			     const LatticeFermion& chi, const LatticeFermion& psi, 
			     enum PlusMinus isign) const
  {
    QDP_error_exit("Hamber-Wu deriv not correct yet");


    START_CODE();

    // This does both parities
    D.deriv(ds_u, chi, psi, isign);

    // Factor in front of the dslash
    for(int mu = 0; mu < Nd; ++mu)
      ds_u[mu] *= fact2;

    // NOTE: missing derivative of 2 link piece

    getFermBC().zero(ds_u);

    END_CODE();
  }

}; // End Namespace Chroma
