// $Id: prec_asqtad_fermact_s.cc,v 1.1 2003-12-10 14:24:08 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */
// NEW $Id: asqtad_fermact_s.cc 2003/11/12 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
//#include "actions/ferm/fermacts/asqtad_fermact_s.h"
//#include "actions/ferm/linop/lmdagm_s.h"

#include "action/ferm/linop/prec_asq_mdagm_s.h"
#include "actions/ferm/linop/asqtad_linop_s.h"
#include "actions/ferm/fermacts/asqtad_fermact_s.h"


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
EvenOddPrecAsqtadFermAct::linOp(const AsqtadConnectState& state) const
{
  return new AsqtadLinOp(state.getFatLinks(), state.getTripleLinks(), Mass);
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
EvenOddPrecAsqtadFermAct::lMdagM(const AsqtadConnectState& state) const
{
  return new AsqtadMdagM(state.getFatLinks(), state.getTripleLinks(), Mass);
}

