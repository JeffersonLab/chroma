// $Id: prec_asqtad_fermact_s.cc,v 1.3 2003-12-11 17:11:17 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */
// NEW $Id: asqtad_fermact_s.cc 2003/11/12 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
//#include "actions/ferm/fermacts/asqtad_fermact_s.h"
//#include "actions/ferm/linop/lmdagm_s.h"

#include "actions/ferm/linop/prec_asq_mdagm_s.h"
#include "actions/ferm/linop/prec_asqtad_linop_s.h"
#include "actions/ferm/fermacts/prec_asqtad_fermact_s.h"


//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u_fat, u_triple 	 fat7 and triple links    (Read)
 * \u has already had KS phases multiplied in.
 */
const EvenOddPrecLinearOperator<LatticeFermion>* 
EvenOddPrecAsqtadFermAct::linOp(const ConnectState& state_) const
{

  // Why in fact are we casting to the base class on both sides of
  // this assignment ? The answer is so that we can use a proxy.
  // Both the Proxy and the ConnectState inherit from the BaseClass
  // and can be cast to and from the base class. However the Proxy
  // and the connect state cannot be directly cast into each other.
  // Which is why we have a virtual base class in the first place.
  //
  // So We cast the ConnectState to an AsqtadConnectStateBase
  // This we can do at our leisure from either AsqtadConnectState
  // OR from the Proxy. We then get access to all the virtual methods
  // in the AsqtadConnectState. Only Restriction: We have to use the
  // get() methods as they are all the base class provides.
  const AsqtadConnectStateBase<LatticeFermion>& state = 
    dynamic_cast<const AsqtadConnectStateBase<LatticeFermion>&>(state_);

  return new EvenOddPrecAsqtadLinOp(state.getFatLinks(), state.getTripleLinks(), Mass);
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
EvenOddPrecAsqtadFermAct::lMdagM(const ConnectState& state_) const
{
  const AsqtadConnectStateBase<LatticeFermion>& state = 
    dynamic_cast<const AsqtadConnectStateBase<LatticeFermion>&>(state_);
  
  return new PrecAsqtadMdagM(state.getFatLinks(), state.getTripleLinks(), Mass);
}

