// $Id: unprec_clover_fermact_w.cc,v 1.3 2003-12-02 15:45:04 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "actions/ferm/fermacts/unprec_clover_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*!
 * \param Mass_        fermion kappa    (Read)
 * \param ClovCoeff_   clover coeff.    (Read)
 * \param u0_          u0    (Read)
 */
void UnprecCloverFermAct::create(const Real& Mass_, const Real& ClovCoeff_, const Real& u0_)
{
  Mass = Mass_;
  ClovCoeff = ClovCoeff_;
  u0   = u0_;
}

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
UnprecCloverFermAct::linOp(const ConnectState& state) const
{
  return new UnprecCloverLinOp(state.getLinks(),Mass,ClovCoeff,u0);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
UnprecCloverFermAct::lMdagM(const ConnectState& state) const
{
  return new lmdagm<LatticeFermion>(UnprecCloverLinOp(state.getLinks(),Mass,ClovCoeff,u0));
}

