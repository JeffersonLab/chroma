// $Id: eoprec_slrc_fermact_w.cc,v 1.1 2007-10-09 03:06:51 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action (fat-relevant, thin-irrelevant terms)
 *
 * Here, the relevant terms are smeared and the irrelevant terms are not smeared.
 * Code provided by Thomas Kaltenbrunner.
 *
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/linop/eoprec_slrc_linop_w.h"
#include "actions/ferm/fermacts/eoprec_slrc_fermact_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"

#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecSLRCFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EvenOddPrecSLRCFermAct(CreateFermStateEnv::reader(xml_in, path), 
					CloverFermActParams(xml_in, path));
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
    const std::string name = "SLRC";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateFermStateEnv::registerAll();
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
   * \param state	    gauge field     	       (Read)
   */
  EvenOddPrecLogDetLinearOperator<LatticeFermion,
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >* 
  EvenOddPrecSLRCFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecSLRCLinOp(state,param);
  }

}

