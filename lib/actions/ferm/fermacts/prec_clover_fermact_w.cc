// $Id: prec_clover_fermact_w.cc,v 1.2 2003-12-02 15:45:04 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_clover_linop_w.h"
#include "actions/ferm/fermacts/prec_clover_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*!
 * \param Mass_        fermion kappa    (Read)
 * \param ClovCoeff_   clover coeff.    (Read)
 * \param u0_          u0    (Read)
 */
void EvenOddPrecCloverFermAct::create(const Real& Mass_, const Real& ClovCoeff_, const Real& u0_)
{
  Mass = Mass_;
  ClovCoeff = ClovCoeff_;
  u0   = u0_;
}

//! Produce a linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>
EvenOddPrecCloverFermAct::linOp(const ConnectState& state) const
{
  return LinearOperator<LatticeFermion>(new EvenOddPrecCloverLinOp(state.getLinks(),Mass,ClovCoeff,u0));
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>
EvenOddPrecCloverFermAct::lMdagM(const ConnectState& state) const
{
  return LinearOperator<LatticeFermion>(EvenOddPrecCloverLinOp(state.getLinks(),Mass,ClovCoeff,u0));
}


//! Computes the derivative of the fermionic action respect to the link field
/*!
 *         |  dS      dS_f
 * ds_u -- | ----   + -----   ( Write )
 *         |  dU       dU
 *
 * psi -- [1./(M_dag*M)]*chi_  ( read ) 
 *
 * \param ds_u     result      ( Write )
 * \param state    gauge field ( Read )
 * \param psi      solution to linear system ( Read )
 */

void
EvenOddPrecCloverFermAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
			       const ConnectState& state,
			       const LatticeFermion& psi) const
{
  START_CODE("EvenOddPrecCloverFermAct::dsdu");
  
  QDPIO::cerr << "EvenOddPrecCloverFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("EvenOddPrecCloverFermAct::dsdu");
}
