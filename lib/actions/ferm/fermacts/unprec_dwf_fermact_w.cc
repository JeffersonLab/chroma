// $Id: unprec_dwf_fermact_w.cc,v 1.6 2004-01-07 13:50:07 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "fermacts.h"
#include "actions/ferm/linop/linop.h"

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::linOp(Handle<const ConnectState> state) const
{
  return new UnprecDWLinOp(state->getLinks(),WilsonMass,m_q);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<LatticeDWFermion>(linOp(state));
}

//! Produce a linear operator for this action but with quark mass 1
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::linOpPV(Handle<const ConnectState> state) const
{
  return new UnprecDWLinOp(state->getLinks(),WilsonMass,1.0);  // fixed to quark mass 1
}

