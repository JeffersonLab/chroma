// $Id: asqtad_fermact_s.cc,v 3.2 2006-11-17 02:17:31 edwards Exp $
/*! \file
 *  \brief Asqtad staggered fermion action
 */
// NEW $Id: asqtad_fermact_s.cc 2003/11/12 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
//#include "actions/ferm/fermacts/asqtad_fermact_s.h"
//#include "actions/ferm/linop/lmdagm_s.h"

#include "actions/ferm/fermacts/fermact_factory_s.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_s.h"

#include "actions/ferm/linop/asqtad_mdagm_s.h"
#include "actions/ferm/linop/asqtad_linop_s.h"
#include "actions/ferm/fermacts/asqtad_fermact_s.h"
#include "util/gauge/stag_phases_s.h"

namespace Chroma 
{ 

  //! Hooks to register the class with the fermact factory
  namespace AsqtadFermActEnv
  {
    //! Callback function
    StaggeredTypeFermAct<LatticeStaggeredFermion,
			 multi1d<LatticeColorMatrix>,
			 multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
									const std::string& path)
    {
      return new AsqtadFermAct(StaggeredCreateFermStateEnv::reader(xml_in, path), 
			       AsqtadFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeStaggeredFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "ASQTAD";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheStagFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheStagTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
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
  EvenOddLinearOperator<LatticeStaggeredFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* 
  AsqtadFermAct::linOp(Handle< FermState<T,P,Q> > state) const
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
    return new AsqtadLinOp(state.cast<AsqtadConnectStateBase>(), param.Mass);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the checkerboarded lattice
   *
   * \param u_fat, u_triple 	 fat7 and triple links	       (Read)
   */
  DiffLinearOperator<LatticeStaggeredFermion, 
		     multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> >* 
  AsqtadFermAct::lMdagM(Handle< FermState<T,P,Q> > state) const
  {
    return new AsqtadMdagM(state.cast<AsqtadConnectStateBase>(), param.Mass);
  }


  //! Create a state -- this multiplies in the K-S phases computes the fat and triple links etc
  AsqtadConnectStateBase*
  AsqtadFermAct::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    multi1d<LatticeColorMatrix> u_with_phases(Nd);
    multi1d<LatticeColorMatrix> u_fat(Nd);
    multi1d<LatticeColorMatrix> u_triple(Nd);

    // First put in the BC
    u_with_phases = u_;
    getFermBC().modify(u_with_phases);

    // Now multiply in phases
    //
    // alpha comes from the StagPhases:: namespace
    for(int i = 0; i < Nd; i++) { 
      u_with_phases[i] *= StagPhases::alpha(i);
    }

    // Make Fat7 and triple links
    Fat7_Links(u_with_phases, u_fat, param.u0);
    Triple_Links(u_with_phases, u_triple, param.u0);

    return new AsqtadConnectState(cfs->getFermBC(), u_with_phases, u_fat, u_triple);
  }

} // End Namespace Chroma

