// $Id: unprec_nef_fermact_array_w.cc,v 1.2 2004-08-24 20:56:13 kostas Exp $
/*! \file
 *  \brief Unpreconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >* 
UnprecNEFFermActArray::linOp(Handle<const ConnectState> state) const
{
  return new UnprecNEFDWLinOpArray(state->getLinks(),WilsonMass,b5,c5,m_q,N5);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >* 
UnprecNEFFermActArray::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
}

//! Produce a linear operator for this action but with quark mass 1
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >* 
UnprecNEFFermActArray::linOpPV(Handle<const ConnectState> state) const
{
  return new UnprecNEFDWLinOpArray(state->getLinks(),WilsonMass,b5,c5,1.0,N5);  // fixed to quark mass 1
}

