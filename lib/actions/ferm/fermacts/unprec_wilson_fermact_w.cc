// $Id: unprec_wilson_fermact_w.cc,v 1.5 2003-10-20 20:31:50 edwards Exp $
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
 * \param _Kappa   fermion kappa    (Read)
 */
void UnprecWilsonFermAct::create(const Real& _Kappa)
{
  Kappa = _Kappa;
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
UnprecWilsonFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u) const
{
  LinearOperator<LatticeFermion>* mdagm = new lmdagm<LatticeFermion>(UnprecWilsonLinOp(u,Kappa));
  return mdagm;
}

