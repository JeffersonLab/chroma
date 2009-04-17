// $Id: eoprec_wilson_coarse_fine_fermact_w.cc,v 3.2 2009-04-17 02:05:30 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action supporting 2+2 anisotropy
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_wilson_coarse_fine_fermact_w.h"
#include "actions/ferm/linop/eoprec_wilson_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"

#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecWilsonCoarseFineFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion, 
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EvenOddPrecWilsonCoarseFineFermAct(CreateFermStateEnv::reader(xml_in, path), 
						    WilsonCoarseFineFermActParams(xml_in, path));
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
    const std::string name("WILSON_COARSE_FINE");

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


  // Init params
  void EvenOddPrecWilsonCoarseFineFermAct::init(const WilsonCoarseFineFermActParams& p)
  {
    param.Mass = p.Mass;
    param.coeffs.resize(Nd);
    param.coeffs = 1.0;

    Real inv_gamma_f = Real(1) / p.gamma_f;

    for(int mu=0; mu < Nd; ++mu)
    {
      if (p.coarse_dirs[mu])
	param.coeffs[mu] = inv_gamma_f;
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
  EvenOddPrecWilsonCoarseFineFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecWilsonLinOp(state,param.Mass,param.coeffs);
  }


  //! Return a linear operator solver for this action to solve M*psi=chi 
  LinOpSystemSolver<LatticeFermion>* 
  EvenOddPrecWilsonCoarseFineFermAct::invLinOp(Handle< FermState<T,P,Q> > state,
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
  EvenOddPrecWilsonCoarseFineFermAct::invMdagM(Handle< FermState<T,P,Q> > state,
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
