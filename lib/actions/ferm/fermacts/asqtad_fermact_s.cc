// $Id: asqtad_fermact_s.cc,v 2.1 2006-01-12 05:45:16 edwards Exp $
/*! \file
 *  \brief Asqtad staggered fermion action
 */
// NEW $Id: asqtad_fermact_s.cc 2003/11/12 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
//#include "actions/ferm/fermacts/asqtad_fermact_s.h"
//#include "actions/ferm/linop/lmdagm_s.h"

#include "actions/ferm/linop/asqtad_mdagm_s.h"
#include "actions/ferm/linop/asqtad_linop_s.h"
#include "actions/ferm/fermacts/asqtad_state.h"
#include "actions/ferm/fermacts/asqtad_fermact_s.h"
#include "util/gauge/stag_phases_s.h"

namespace Chroma 
{ 

//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u_fat, u_triple 	 fat7 and triple links    (Read)
 * \u has already had KS phases multiplied in.
 */
const EvenOddLinearOperator< LatticeStaggeredFermion, multi1d<LatticeColorMatrix> >* 
AsqtadFermAct::linOp(Handle<const ConnectState> state_) const
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
  const AsqtadConnectStateBase<LatticeStaggeredFermion>& state = 
    dynamic_cast<const AsqtadConnectStateBase<LatticeStaggeredFermion>&>(*state_);

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

const DiffLinearOperator<LatticeStaggeredFermion, multi1d<LatticeColorMatrix> >* 
AsqtadFermAct::lMdagM(Handle<const ConnectState> state_) const
{
  const AsqtadConnectStateBase<LatticeStaggeredFermion>& state = 
    dynamic_cast<const AsqtadConnectStateBase<LatticeStaggeredFermion>&>(*state_);
  
  return new AsqtadMdagM(state.getFatLinks(), state.getTripleLinks(), Mass);
}

//! Create a state -- this multiplies in the K-S phases computes the fat and triple links etc


const AsqtadConnectStateBase<LatticeStaggeredFermion>*
AsqtadFermAct::createState(const multi1d<LatticeColorMatrix>& u_) const
{
  multi1d<LatticeColorMatrix> u_with_phases(Nd);
  multi1d<LatticeColorMatrix> u_fat(Nd);
  multi1d<LatticeColorMatrix> u_triple(Nd);

  // First put in the BC
  u_with_phases = u_;
  getFermBC().modifyU(u_with_phases);

  // Now multiply in phases
  //
  // alpha comes from the StagPhases:: namespace
  for(int i = 0; i < Nd; i++) { 
    u_with_phases[i] *= StagPhases::alpha(i);
  }

  // Make Fat7 and triple links
  Fat7_Links(u_with_phases, u_fat, u0);
  Triple_Links(u_with_phases, u_triple, u0);

  return new AsqtadConnectState<LatticeStaggeredFermion>(u_with_phases, u_fat, u_triple);
}

}; // End Namespace Chroma

