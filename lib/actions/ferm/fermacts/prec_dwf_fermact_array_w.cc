// $Id: prec_dwf_fermact_array_w.cc,v 1.5 2004-01-07 13:50:07 bjoo Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"
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
EvenOddPrecDWFermActArray::linOp(Handle<const ConnectState> state) const
{
  return new EvenOddPrecDWLinOpArray(state->getLinks(),WilsonMass,m_q,N5);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecDWFermActArray::lMdagM(Handle<const ConnectState> state) const
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
EvenOddPrecDWFermActArray::linOpPV(Handle<const ConnectState> state) const
{
  // For the PV operator, use the **unpreconditioned** one
  // fixed to quark mass 1
  return new UnprecDWLinOpArray(state->getLinks(),WilsonMass,1.0,N5);
}

