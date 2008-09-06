// $Id: klein_gordon_fermact_s.cc,v 3.2 2008-09-06 18:36:25 bjoo Exp $
/*! \file
 *  \brief Klein-Gordon boson action masquerading action as a staggered action
 */

#include "actions/ferm/fermacts/fermact_factory_s.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_s.h"

#include "actions/ferm/linop/klein_gordon_linop_s.h"
#include "actions/ferm/fermacts/klein_gordon_fermact_s.h"

namespace Chroma 
{ 
  // Default constructor
  KleinGordonFermActParams::KleinGordonFermActParams()
  {
    Mass = 0.0;
  }

  // Read parameters
  KleinGordonFermActParams::KleinGordonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "Mass", Mass);

    //  Read optional anisoParam.
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, KleinGordonFermActParams& param)
  {
    KleinGordonFermActParams tmp(xml, path);
    param = tmp;
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const KleinGordonFermActParams& param)
  {
    push(xml, path);

    write(xml, "Mass", param.Mass);
    if (param.anisoParam.anisoP)
      write(xml, "AnisoParam", param.anisoParam);

    pop(xml);
  }


  //! Hooks to register the class with the fermact factory
  namespace KleinGordonFermActEnv
  {
    //! Callback function
    StaggeredTypeFermAct<LatticeStaggeredFermion,
			 multi1d<LatticeColorMatrix>,
			 multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
									const std::string& path)
    {
      return new KleinGordonFermAct(StaggeredCreateFermStateEnv::reader(xml_in, path), 
				    KleinGordonFermActParams(xml_in, path));
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
    const std::string name = "KLEIN_GORDON";

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


  //! Create a state -- this multiplies in the K-S phases computes the fat and triple links etc
  FermState<LatticeStaggeredFermion, 
	    multi1d<LatticeColorMatrix>,
	    multi1d<LatticeColorMatrix> >* 
  KleinGordonFermAct::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    return getCreateState()(u_);
  }


  // Produce a linear operator for this action
  UnprecLinearOperator<LatticeStaggeredFermion,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >* 
  KleinGordonFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new KleinGordonLinOp(state, param.Mass, param.anisoParam);
  }


  // Produce a M^dag.M linear operator for this action
  DiffLinearOperator<LatticeStaggeredFermion, 
		     multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix> >* 
  KleinGordonFermAct::lMdagM(Handle< FermState<T,P,Q> > state) const
  {
    return new DiffMdagMLinOp<T,P,Q>(this->linOp(state));
  }


#if 0
  // Already supplied in  chroma/lib/actions/ferm/qprop/quarkprop_s.cc

  // Return a linear operator solver for this action to solve MdagM*psi=chi 
  MdagMSystemSolver<LatticeStaggeredFermion>* 
  KleinGordonFermAct::invMdagM(Handle< FermState<T,P,Q> > state,
			       const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMStagFermSystemSolverFactory::Instance().createObject(invParam.id,
									paramtop,
									invParam.path,
									linOp(state));
  }


  // Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
  MdagMMultiSystemSolver<LatticeStaggeredFermion>* 
  KleinGordonFermAct::mInvMdagM(Handle< FermState<T,P,Q> > state,
				const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMStagFermMultiSystemSolverFactory::Instance().createObject(invParam.id,
									     paramtop,
									     invParam.path,
									     lMdagM(state));
  }

  // Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
  MdagMMultiSystemSolverAccumulate<LatticeStaggeredFermion>* 
  KleinGordonFermAct::mInvMdagMAcc(Handle< FermState<T,P,Q> > state,
				const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMStagFermMultiSystemSolverAccumulateFactory::Instance().createObject(invParam.id,
									     paramtop,
									     invParam.path,
									     lMdagM(state));
  }


  // Return quark prop solver, solution of unpreconditioned system
//  SystemSolver<T>* 
//  KleinGordonFermAct::qprop(Handle< FermState<T,P,Q> > state,
//			    const GroupXML_t& invParam) const;
#endif


} // End Namespace Chroma

