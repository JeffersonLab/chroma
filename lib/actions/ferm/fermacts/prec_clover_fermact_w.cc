// $Id: prec_clover_fermact_w.cc,v 1.3 2004-01-02 03:19:40 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_clover_linop_w.h"
#include "actions/ferm/fermacts/prec_clover_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Produce a linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state	    gauge field     	       (Read)
 */
const EvenOddPrecLinearOperator<LatticeFermion>* 
EvenOddPrecCloverFermAct::linOp(Handle<const ConnectState> state) const
{
  return new EvenOddPrecCloverLinOp(state->getLinks(),Mass,ClovCoeff,u0);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
EvenOddPrecCloverFermAct::lMdagM(Handle<const ConnectState> state) const
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
EvenOddPrecCloverFermAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
			       Handle<const ConnectState> state,
			       const LatticeFermion& psi) const
{
  START_CODE("EvenOddPrecCloverFermAct::dsdu");
  
  QDPIO::cerr << "EvenOddPrecCloverFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("EvenOddPrecCloverFermAct::dsdu");
}
