// $Id: unprec_wilson_fermact_w.cc,v 1.19 2004-08-05 15:59:03 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"


//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>*
UnprecWilsonFermAct::linOp(Handle<const ConnectState> state) const
{
  return new UnprecWilsonLinOp(state->getLinks(), Mass); 
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>*
UnprecWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<LatticeFermion>(linOp(state));
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
  for(int mu = 0; mu < Nd; mu++) { 

    // f_tmp = (1 + gamma_mu) Y 
    f_tmp = Gamma(1<<mu)*Y;
    f_tmp += Y;

    //   trace_spin ( X^{dag}_x ( 1 + gamma_mu ) Y_x+mu )
    // = trace_spin(  Y_x+mu X^dag_x (1 + gamma_mu )  )
    u_tmp = trace(adj(psi)*shift(f_tmp, FORWARD, mu));

    // f_tmp = -(1 -gamma_mu) X
    f_tmp = Gamma(1<<mu)*psi;
    f_tmp -= psi;

    //  -trace_spin( -Y^{dag}_x( 1 - gamma_mu) X_x+mu )
    //= +trace_spin( X_{x+mu} Y^{dag}_x ( 1 - gamma_mu ) 
    u_tmp -= trace(adj(Y)*shift(f_tmp, FORWARD, mu));
    
    // accumulate with prefactor
    ds_u[mu] += prefactor*u[mu]*u_tmp;
  }

     
    
  END_CODE();
}
