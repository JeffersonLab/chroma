// $Id: unprec_wilson_fermact_w.cc,v 1.12 2003-12-15 17:52:51 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*!
 * \param Mass_   fermion kappa    (Read)
 */
void UnprecWilsonFermAct::create(const Real& Mass_)
{
  Mass = Mass_;
//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}

//! Produce a linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>*
UnprecWilsonFermAct::linOp(const ConnectState& state) const
{
  const UnprecWilsonLinOp* D = new UnprecWilsonLinOp(state.getLinks(), Mass);
  return D;

}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the entire lattice
 *
 * \param state    gauge field     	       (Read)
 */
const LinearOperator<LatticeFermion>*
UnprecWilsonFermAct::lMdagM(const ConnectState& state) const
{

 const UnprecWilsonLinOp* D = new UnprecWilsonLinOp(state.getLinks(), Mass);
  return new lmdagm<LatticeFermion>(*D);
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
UnprecWilsonFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
			  const ConnectState& state,
			  const LatticeFermion& psi) const
{
  START_CODE("UnprecWilsonFermAct::dsdu");
  
  QDPIO::cerr << "UnprecWilsonFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("UnprecWilsonFermAct::dsdu");
}
