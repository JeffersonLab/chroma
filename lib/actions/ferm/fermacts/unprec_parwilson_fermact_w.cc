// $Id: unprec_parwilson_fermact_w.cc,v 1.1 2004-01-12 04:48:00 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"
#include "actions/ferm/linop/unprec_parwilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"


//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>*
UnprecParWilsonFermAct::linOp(Handle<const ConnectState> state) const
{
  return new UnprecParWilsonLinOp(state->getLinks(),Mass,H); 
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>*
UnprecParWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
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
UnprecParWilsonFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
			     Handle<const ConnectState> state,
			     const LatticeFermion& psi) const
{
  START_CODE("UnprecParWilsonFermAct::dsdu");
  
  QDPIO::cerr << "UnprecParWilsonFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("UnprecParWilsonFermAct::dsdu");
}
