// $Id: prec_slic_fermact_w.cc,v 3.2 2006-09-19 16:04:22 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_slic_linop_w.h"
#include "actions/ferm/fermacts/prec_slic_fermact_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"

#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecSLICFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
#if 1
      return new EvenOddPrecSLICFermAct(CreateFermStateEnv::reader(xml_in, path), 
					  CloverFermActParams(xml_in, path));
#else 
      Handle< CreateFermState<LatticeFermion,
	multi1d<LatticeColorMatrix>,
	multi1d<LatticeColorMatrix> > > cfs(
	  new CreateStoutFermState(WilsonTypeFermBCEnv::reader(xml_in, path + "/FermState"),
				   StoutFermStateParams(xml_in, path + "/FermState")));
      
      return new EvenOddPrecSLICFermAct(cfs,
					  CloverFermActParams(xml_in, path));
#endif
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
    const std::string name = "SLIC";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= CreateFermStateEnv::registered;

      foo &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
      foo &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Produce a linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state	    gauge field     	       (Read)
   */
  EvenOddPrecLogDetLinearOperator<LatticeFermion,
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >* 
  EvenOddPrecSLICFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecSLICLinOp(state,param);
  }

}

