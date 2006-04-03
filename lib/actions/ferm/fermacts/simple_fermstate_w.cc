// $Id: simple_fermstate_w.cc,v 3.0 2006-04-03 04:58:46 edwards Exp $
/*! \file
 *  \brief Simple ferm state and a creator
 */

#include "actions/ferm/fermacts/simple_fermstate.h"
#include "actions/ferm/fermacts/simple_fermstate_w.h"
#include "actions/ferm/fermacts/ferm_createstate_factory_w.h"
#include "actions/ferm/fermacts/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{

  /*! \ingroup fermacts */
  namespace CreateSimpleFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
							      const std::string& path) 
    {
      return new CreateSimpleFermState<LatticeFermion,
	                               multi1d<LatticeColorMatrix>, 
	                               multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml, 
												 path));
    }

    const std::string name = "SIMPLE_FERM_STATE";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, 
									  createFerm);
      return foo;
    }

    const bool registered = registerAll();
  }

}

