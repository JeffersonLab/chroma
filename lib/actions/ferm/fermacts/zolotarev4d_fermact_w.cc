// $Id: zolotarev4d_fermact_w.cc,v 1.3 2003-12-02 15:45:04 edwards Exp $
/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/zolotarev4d_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param _m_q    quark mass   	       (Read)
 * \param _M 	  operator kernel      (Read)
 */
void Zolotarev4DFermAct::Zolotarev4DFermAct(const Real& _m_q, const UnprecWilsonTypeFermAct& _M) :
  m_q(_m_q), M(_M) 
{
  RsdCGinner = 1.0e-7;  // Hardwired the accuracy

  NEigVal = 0;
}

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermAct::linOp(const ConnectState& state) const
{
  return new UnprecWilsonLinOp(state.getLinks(),Kappa);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * Should use special form when we know we have exact chiral symmetry
 *
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermAct::lMdagM(const ConnectState& state) const
{
  //  *****NOTE***** 
  // Should use special form when we know we have exact chiral symmetry
  return new lmdagm(UnprecWilsonLinOp(state.getLinks(),Kappa));
}

