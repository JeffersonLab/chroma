// $Id: asqtad_fermact_s.cc,v 1.1 2003-12-10 12:38:14 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */
// NEW $Id: asqtad_fermact_s.cc 2003/11/12 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
//#include "actions/ferm/fermacts/asqtad_fermact_s.h"
//#include "actions/ferm/linop/lmdagm_s.h"

#include "prec_asq_mdagm_s.h"
#include "asqtad_linop_s.h"
#include "asqtad_fermact_s.h"
//#include "lmdagm_s.h"

//! Creation routine
/*! \ingroup fermact
 *
 * \param _Mass   fermion mass    (Read)
 */
void AsqtadFermAct::create(const Real& _Mass)
{
  Mass = _Mass;
//    CoeffWilsr_s = (AnisoP) ? Wilsr_s / xiF_0 : 1;
}

//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u_fat, u_triple 	 fat7 and triple links    (Read)
 * \u has already had KS phases multiplied in.
 */
const LinearOperator<LatticeFermion>* 
AsqtadFermAct::linOp(const multi1d<LatticeColorMatrix>& u_fat, const multi1d<LatticeColorMatrix>& u_triple) const
{
  return new AsqtadLinOp(u_fat,u_triple,Mass);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the checkerboarded lattice
 *
 * \param u_fat, u_triple 	 fat7 and triple links	       (Read)
 */
const LinearOperator<LatticeFermion>* 
AsqtadFermAct::lMdagM(const multi1d<LatticeColorMatrix>& u_fat, const 
multi1d<LatticeColorMatrix>& u_triple) const
{
  return new AsqtadMdagM(u_fat, u_triple, Mass);
}

