// $Id: unprec_dwf_fermact_w.cc,v 1.3 2003-11-09 22:33:01 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_dwf_linop_w.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param WilsonMass   DWF height    (Read)
 * \param m_q          quark mass    (Read)
 */
void UnprecDWFermAct::create(const Real& WilsonMass_, const Real& m_q_)
{
  WilsonMass = WilsonMass_;
  m_q = m_q_;
  a5  = 1.0;
//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}

//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::linOp(const multi1d<LatticeColorMatrix>& u) const
{
  return new UnprecDWLinOp(u,WilsonMass,m_q);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<LatticeDWFermion>* mdagm = new lmdagm<LatticeDWFermion>(UnprecDWLinOp(u,WilsonMass,m_q));
  return mdagm;
}

//! Produce a linear operator for this action but with quark mass 1
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeDWFermion>* 
UnprecDWFermAct::linOpPV(const multi1d<LatticeColorMatrix>& u) const
{
  return new UnprecDWLinOp(u,WilsonMass,1.0);  // fixed to quark mass 1
}

