// $Id: prec_ovdwf_fermact_array_w.cc,v 1.1 2004-02-13 20:56:43 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned Overlap-DWF (Borici) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_ovdwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovdwf_linop_array_w.h"
#include "actions/ferm/linop/prec_ovdwf_linop_array_w.h"
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
EvenOddPrecOvDWFermActArray::linOp(Handle<const ConnectState> state) const
{
  return new EvenOddPrecOvDWLinOpArray(state->getLinks(),WilsonMass,m_q,N5);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecOvDWFermActArray::lMdagM(Handle<const ConnectState> state) const
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
EvenOddPrecOvDWFermActArray::linOpPV(Handle<const ConnectState> state) const
{
  // For the PV operator, use the **unpreconditioned** one
  // fixed to quark mass 1
  return new UnprecDWLinOpArray(state->getLinks(),WilsonMass,1.0,N5);
}

