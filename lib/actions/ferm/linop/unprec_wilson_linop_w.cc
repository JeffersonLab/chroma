// $Id: unprec_wilson_linop_w.cc,v 1.12 2004-12-17 17:45:04 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma 
{ 

  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_ 	    gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   */
  void UnprecWilsonLinOp::create(const multi1d<LatticeColorMatrix>& u_, const Real& Mass_)
  {
    Mass = Mass_;
    u = u_;
    D.create(u);

    //    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
  }


  //! Apply unpreconditioned Wilson fermion linear operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field     	       (Read)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecWilsonLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				      enum PlusMinus isign) const
  {
    START_CODE();

    //
    //  Chi   =  (Nd+Mass)*Psi  -  (1/2) * D' Psi
    //
    LatticeFermion tmp;
    Real fact1 = Nd + Mass;
    Real fact2 = -0.5;

    // D is a Dslash - must apply to both CB-s
    D(tmp, psi, isign);

    chi = fact1*psi + fact2*tmp;
  
    END_CODE();
  }



  //! Derivative of unpreconditioned Wilson dM/dU
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of chi vector                  (Read)
   *
   * \return Computes   chi^dag * \dot(D} * psi  
   */
  void 
  UnprecWilsonLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
  {
    START_CODE();

    // This does both parities
    D.deriv(ds_u, chi, psi, isign);

    // Factor from the -1/2 in front of the dslash
    for(int mu = 0; mu < Nd; ++mu)
      ds_u[mu] *= Real(-0.5);

    END_CODE();
  }


#if 0
  //! Derivative of unpreconditioned Wilson dM/dU
  /*!
   * This subroutine applies the derivative of the Unprecontitioned
   * Wilson Operator with respect to the gauge fields to a vector 
   * Psi,
   * 
   *  isign = PLUS: 
   *
   *  Chi = dM/dU 
   *      = -(1/2) [  ( 1 - gamma_mu ) psi(x+mu) ]
   *
   *  isign = MINUS:
   * 
   *  Chi = dM^dagger/dU 
   *      = -(1/2) [  ( 1 + gamma_mu ) psi(x+mu) ]
   *
   * Amusingly enough this is a linearop that does not depend 
   * on the gauge fields....
   *
   * The mu is done at creation...
   *
   */
  void 
  UnprecWilsonLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    LatticeFermion tmp;
    LatticeFermion gamma_mu_psi;
    
    for(int mu = 0; mu < Nd; ++mu)
    {
      // Get gamma_mu * psi
      gamma_mu_psi = Gamma(1 << mu)*psi;

      // Construct the right derivative term...
      if( isign == PLUS ) 
      {
	// Undaggered:
	// ( 1 - gamma_mu ) psi
	tmp = psi - gamma_mu_psi;
      }
      else 
      { 
	// Daggered:
	// ( 1 + gamma_mu) psi
	tmp = psi + gamma_mu_psi;
      }
    
      // Do the shift here for delta_y, x+mu and form derivative
      ds_u[mu] = traceSpin(outerProduct(-Real(0.5)*shift(tmp, FORWARD, mu),chi));
    }
    
    END_CODE();
  }


  //! Computes the derivative of the fermionic action respect to the link field
  /*!
   *         |  dS      dS_f
   * ds_u -- | ----   + -----   ( Write )
   *         |  dU       dU
   *
   * psi -- [1./(M_dag*M)]*chi_  ( read ) 
   *
   * \param ds_u     result      ( Write )
   * \param state    gauge field ( Read )
   * \param psi      solution to linear system ( Read )
   */

  void
  UnprecWilsonFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
			    Handle<const ConnectState> state,
			    const LatticeFermion& psi) const
  {
    START_CODE();

    // The phi <=> X so I define the Y field as MX 
  
    // Get at the U matrices
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
  
    // Get a linear operator
    Handle<const LinearOperator<LatticeFermion> > M(linOp(state));

    // Compute MY
    LatticeFermion Y;
    (*M)(Y, psi, PLUS);

    // Usually this is Kappa. In our normalisation it is 0.5 
    // I am adding in a factor of -1 to be consistent with the sign
    // convention for the preconditioned one. (We can always take this out
    // later
    Real prefactor=-Real(0.5);

    // Two temporaries
    LatticeFermion f_tmp;
    LatticeColorMatrix u_tmp;
    for(int mu = 0; mu < Nd; mu++)
    { 
      // f_tmp = (1 + gamma_mu) Y 
      f_tmp = Gamma(1<<mu)*Y;
      f_tmp += Y;

      //   trace_spin ( ( 1 + gamma_mu ) Y_x+mu X^{dag}_x )
//      u_tmp = traceSpin(outerProduct(shift(f_tmp, FORWARD, mu),psi));
      LatticeFermion foo = shift(f_tmp, FORWARD, mu);
      u_tmp = traceSpin(outerProduct(foo,psi));

      // f_tmp = -(1 -gamma_mu) X
      f_tmp = Gamma(1<<mu)*psi;
      f_tmp -= psi;

      //  +trace_spin( ( 1 - gamma_mu) X_x+mu Y^{dag}_x)
//      u_tmp -= traceSpin(outerProduct(shift(f_tmp, FORWARD, mu),Y));
      foo = shift(f_tmp, FORWARD, mu);
      u_tmp -= traceSpin(outerProduct(foo,Y));
    
      // accumulate with prefactor
      ds_u[mu] = prefactor*( u[mu]*u_tmp );
    }
    
    END_CODE();
  }

  void
  UnprecWilsonFermAct::dsdu2(multi1d<LatticeColorMatrix>& ds_u,
			    Handle<const ConnectState> state,
			    const LatticeFermion& psi) const
  {
    START_CODE();

    // The phi <=> X so I define the Y field as MX 
  
    // Get at the U matrices
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
  
    // Get a linear operator
    Handle<const LinearOperator<LatticeFermion> > M(linOp(state));

    // Compute MY
    LatticeFermion Y;
    (*M)(Y, psi, PLUS);

    Real prefactor=-Real(0.5);

    // Derivative Matrix
    for(int mu=0; mu < Nd; mu++) { 

      // Create a linop for dM/dU_mu
      Handle<const LinearOperator<LatticeFermion> > dMdU(ldMdU(state,mu));

      LatticeColorMatrix u_tmp;
      LatticeFermion tmp;

      // tmp = (dM^{dagger}/dU_mu) Y
      (*dMdU)(tmp, Y, MINUS);

      // u_tmp += Tr_spin (dM^{dagger}/dU_mu) Y X^{dagger} 
      //        = traceSpin(outerProduct(tmp, X))
      u_tmp = traceSpin(outerProduct(tmp, psi));

    
      // tmp = (dM/dU_mu) X = (dM/dU_w) phi 
      (*dMdU)(tmp, psi, PLUS);

      // u_tmp = Tr_spin (dM/dU_mu) X Y^{dagger} = traceSpin(outerProduct(tmp, Y))
      u_tmp += traceSpin(outerProduct(tmp, Y));

    
      // Multiply in u_[mu]
      ds_u[mu] = u[mu]*u_tmp;

      
      

    }
    
    END_CODE();
  }
#endif

}; // End Namespace Chroma
