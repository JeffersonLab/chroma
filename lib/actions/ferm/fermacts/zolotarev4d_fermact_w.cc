// $Id: zolotarev4d_fermact_w.cc,v 1.4 2003-12-02 22:35:26 edwards Exp $
/*! \file
 *  \brief 4D Zolotarev variant of Overlap-Dirac operator
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/zolotarev4d_fermact_w.h"
#include "actions/ferm/linop/zolotarev4d_linop_w.h"

// Replace this with special overlap M^dag*M version
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*! */
void 
Zolotarev4DFermAct::init()
{
}

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field state  	       (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermAct::linOp(const EVConnectState<LatticeFermion>& state) const
{
  return new Zolotarev4DLinOp(state,M,m_q);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * Should use special form when we know we have exact chiral symmetry
 *
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field state   	       (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermAct::lMdagM(const EVConnectState<LatticeFermion>& state) const
{
  //  *****NOTE***** 
  // Should use special form when we know we have exact chiral symmetry
  return new lmdagm<LatticeFermion>(Zolotarev4DLinOp(state,M,m_q));
}

