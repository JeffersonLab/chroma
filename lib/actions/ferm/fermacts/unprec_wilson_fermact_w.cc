// $Id: unprec_wilson_fermact_w.cc,v 1.9 2003-11-22 21:32:10 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/lmdagm_w.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param Mass_   fermion kappa    (Read)
 */
void UnprecWilsonFermAct::create(const Real& Mass_)
{
  Mass = Mass_;
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
const LinearOperator<LatticeFermion>* 
UnprecWilsonFermAct::linOp(const multi1d<LatticeColorMatrix>& u) const
{
  return new UnprecWilsonLinOp(u,Mass);
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
UnprecWilsonFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<LatticeFermion>* mdagm = new lmdagm<LatticeFermion>(UnprecWilsonLinOp(u,Mass));
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
UnprecWilsonFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
			  const multi1d<LatticeColorMatrix>& u, 
			  const LatticeFermion& psi) const
{
  START_CODE("UnprecWilsonFermAct::dsdu");
  
  QDPIO::cerr << "UnprecWilsonFermAct::dsdu not implemented" << endl;
  QDP_abort(1);

  END_CODE("UnprecWilsonFermAct::dsdu");
}
