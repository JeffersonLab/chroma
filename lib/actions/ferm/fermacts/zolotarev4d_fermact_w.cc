// $Id: zolotarev4d_fermact_w.cc,v 1.2 2003-10-20 20:31:50 edwards Exp $
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
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermAct::linOp(const multi1d<LatticeColorMatrix>& u) const
{
  return new UnprecWilsonLinOp(u,Kappa);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
Zolotarev4DFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<LatticeFermion>* mdagm = new lmdagm(UnprecWilsonLinOp(u,Kappa));
  return mdagm;
}

