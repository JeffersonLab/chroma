// $Id: unprec_clover_fermact_w.cc,v 1.1 2003-11-22 21:33:24 edwards Exp $
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
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
UnprecCloverFermAct::linOp(const multi1d<LatticeColorMatrix>& u) const
{
  return new UnprecCloverLinOp(u,Mass,ClovCoeff,u0);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param u 	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>* 
UnprecCloverFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<LatticeFermion>* mdagm = 
    new lmdagm<LatticeFermion>(UnprecCloverLinOp(u,Mass,ClovCoeff,u0));
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
UnprecCloverFermAct::dsdu(multi1d<LatticeColorMatrix>& ds_u,
			  const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& psi) const
{
  START_CODE("UnprecCloverFermAct::dsdu");
  
  QDPIO::cerr << "UnprecCloverFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("UnprecCloverFermAct::dsdu");
}
