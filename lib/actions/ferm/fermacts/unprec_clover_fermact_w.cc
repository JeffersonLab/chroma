// $Id: unprec_clover_fermact_w.cc,v 1.4 2004-01-02 03:19:40 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "actions/ferm/fermacts/unprec_clover_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
UnprecCloverFermAct::linOp(Handle<const ConnectState> state) const
{
  return new UnprecCloverLinOp(state->getLinks(),Mass,ClovCoeff,u0);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
UnprecCloverFermAct::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<LatticeFermion>(linOp(state));
}

