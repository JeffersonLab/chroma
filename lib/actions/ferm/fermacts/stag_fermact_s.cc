// $Id: stag_fermact_s.cc,v 3.0 2006-04-03 04:58:46 edwards Exp $
/*! \file
 *  \brief Staggered fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/stag_fermact_s.h"
#include "actions/ferm/linop/stag_linop_s.h"
#include "actions/ferm/linop/stag_mdagm_linop_s.h"

namespace Chroma 
{ 
  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  EvenOddLinearOperator<T,P,Q>*
  StagFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new StagLinOp(state, Mass); 
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state    gauge field     	       (Read)
   */
  LinearOperator<LatticeStaggeredFermion>*
  StagFermAct::lMdagM(Handle< FermState<T,P,Q> > state) const
  {
    return new StagMdagMLinOp(state, Mass);
  }


#if 0

#error "NEEDS MORE CONVERSION"

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
  StagFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
		    Handle<const FermState<T,P,Q> > state,
		    const LatticeStaggeredFermion& psi) const
  {
    LatticeStaggeredFermion u_psi;
    LatticeStaggeredFermion rho;
    LatticeStaggeredFermion u_rho;

    LatticeStaggeredFermion tmp_1;
    Real dummy;
    int mu;
    int cb;

    START_CODE();

    /* rho = Dslash(1<-0) * psi */
    dslash (u, psi, rho, PLUS, 0);

    /* rho = (KappaMD^2)*rho = (KappaMC^2)*Dslash*psi */
    dummy = KappaMD*KappaMD;
    rho = rho * dummy;

    for(mu = 0;mu  < ( Nd); ++mu )
    {
      cb = 0;

      /* tmp_1(x) = rho(x+mu) */
      tmp_1[rb[cb]] = shift(rho, FORWARD, mu);

      /* u_rho(x) = u(x,mu)*rho(x+mu) = u(x,mu)*tmp_1(x) */
      /* Note the KS phase factors are already included in the U's! */
      u_rho = u[mu][cb] * tmp_1;
    
      /* ds_u(x,mu) = - u_rho(x) * psi_dag(x,mu) */
      ds_u[mu][cb] -= u_rho * adj(psi);
        
      cb = 1;
    
      /* tmp_1(x) = psi(x+mu) */
      tmp_1[rb[cb]] = shift(psi, FORWARD, mu);
    
      /* u_psi(x) = u(x,mu)*psi(x+mu) = u(x,mu)*tmp_1(x) */
      /* Note the KS phase factors are already included in the U's! */
      u_psi = u[mu][cb] * tmp_1;
    
      /* ds_u(x,mu) = u_psi(x) * rho_dag(x,mu) */
      ds_u[mu][cb] += u_psi * adj(rho);
    }
      
    END_CODE();
  }

} // End Namespace Chroma

#endif
