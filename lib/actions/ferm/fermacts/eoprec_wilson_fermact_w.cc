// $Id: eoprec_wilson_fermact_w.cc,v 3.2 2009-04-17 02:05:30 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_wilson_fermact_w.h"
#include "actions/ferm/linop/eoprec_wilson_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"

#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion, 
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EvenOddPrecWilsonFermAct(CreateFermStateEnv::reader(xml_in, path), 
					  WilsonFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name("WILSON");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
  }


  //! Produce a linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state 	    gauge field     	       (Read)
   */
  EvenOddPrecConstDetLinearOperator<LatticeFermion,
				    multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> >* 
  EvenOddPrecWilsonFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecWilsonLinOp(state,param.Mass,param.anisoParam);
  }


  //! Return a linear operator solver for this action to solve M*psi=chi 
  LinOpSystemSolver<LatticeFermion>* 
  EvenOddPrecWilsonFermAct::invLinOp(Handle< FermState<T,P,Q> > state,
				     const GroupXML_t& invParam) const
  {
    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
	
    return TheLinOpFermSystemSolverFactory::Instance().createObject(invParam.id,
								    paramtop,
								    invParam.path,
								    state,
								    linOp(state));
  }


  //! Return a linear operator solver for this action to solve M^dag.M*psi=chi 
  MdagMSystemSolver<LatticeFermion>* 
  EvenOddPrecWilsonFermAct::invMdagM(Handle< FermState<T,P,Q> > state,
				     const GroupXML_t& invParam) const
  {
    std::istringstream  is(invParam.xml);
    XMLReader  paramtop(is);
	
    return TheMdagMFermSystemSolverFactory::Instance().createObject(invParam.id,
								    paramtop,
								    invParam.path,
								    state,
								    linOp(state));
  }

}
