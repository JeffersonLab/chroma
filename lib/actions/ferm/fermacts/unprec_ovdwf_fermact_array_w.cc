// $Id: unprec_ovdwf_fermact_array_w.cc,v 1.5 2004-01-07 13:50:07 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovdwf_linop_array_w.h"
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
UnprecOvDWFermActArray::linOp(Handle<const ConnectState> state) const
{
  return new UnprecOvDWLinOpArray(state->getLinks(),WilsonMass,m_q,N5);
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
UnprecOvDWFermActArray::lMdagM(Handle<const ConnectState> state) const
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
UnprecOvDWFermActArray::linOpPV(Handle<const ConnectState> state) const
{
  return new UnprecOvDWLinOpArray(state->getLinks(),WilsonMass,1.0,N5);  // fixed to quark mass 1
}

