// $Id: prec_dwf_fermact_array_w.cc,v 1.3 2003-12-02 15:45:04 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param WilsonMass_   DWF height    (Read)
 * \param m_q_          quark mass    (Read)
 * \param N5_           extent of DW flavor space   (Read)
 */
void EvenOddPrecDWFermActArray::create(const Real& WilsonMass_, const Real& m_q_, int N5_)
{
  WilsonMass = WilsonMass_;
  m_q = m_q_;
  N5  = N5_;

  a5  = 1.0;

//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}


//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const EvenOddPrecLinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecDWFermActArray::linOp(const ConnectState& state) const
{
  return new EvenOddPrecDWLinOpArray(state.getLinks(),WilsonMass,m_q,N5);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecDWFermActArray::lMdagM(const ConnectState& state) const
{
  return new lmdagm<multi1d<LatticeFermion> >(EvenOddPrecDWLinOpArray(state.getLinks(),WilsonMass,m_q,N5));
}

//! Produce a linear operator for this action but with quark mass 1
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
EvenOddPrecDWFermActArray::linOpPV(const ConnectState& state) const
{
  // For the PV operator, use the **unpreconditioned** one
  // fixed to quark mass 1
  return new UnprecDWLinOpArray(state.getLinks(),WilsonMass,1.0,N5);
}

