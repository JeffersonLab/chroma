// $Id: prec_nef_fermact_array_w.cc,v 1.2 2004-08-24 20:56:13 kostas Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/prec_nef_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const EvenOddPrecLinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecNEFFermActArray::linOp(Handle<const ConnectState> state) const
{
 return new EvenOddPrecNEFDWLinOpArray(state->getLinks(),WilsonMass,b5,c5,m_q,N5);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecNEFFermActArray::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
}

//! Produce a linear operator for this action but with quark mass 1
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecNEFFermActArray::linOpPV(Handle<const ConnectState> state) const
{
  // For the PV operator, use the **unpreconditioned** one
  // fixed to quark mass 1
  return new UnprecNEFDWLinOpArray(state->getLinks(),WilsonMass,b5,c5,1.0,N5);
}

