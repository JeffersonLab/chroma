// $Id: prec_clover_fermact_w.cc,v 1.1 2003-11-22 21:34:01 edwards Exp $
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
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
EvenOddPrecCloverFermAct::linOp(const multi1d<LatticeColorMatrix>& u) const
{
  return new EvenOddPrecCloverLinOp(u,Mass,ClovCoeff,u0);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd subset
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
EvenOddPrecCloverFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<LatticeFermion>* mdagm = 
    new lmdagm<LatticeFermion>(EvenOddPrecCloverLinOp(u,Mass,ClovCoeff,u0));
  return mdagm;
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
 * \param u        gauge field ( Read )
 * \param psi      solution to linear system ( Read )
 */

void
EvenOddPrecCloverFermAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
			       const multi1d<LatticeColorMatrix>& u, 
			       const LatticeFermion& psi) const
{
  START_CODE("EvenOddPrecCloverFermAct::dsdu");
  
  QDPIO::cerr << "EvenOddPrecCloverFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("EvenOddPrecCloverFermAct::dsdu");
}
