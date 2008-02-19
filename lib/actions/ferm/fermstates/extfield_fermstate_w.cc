// $Id: extfield_fermstate_w.cc,v 1.6 2008-02-19 21:11:29 edwards Exp $
/*! \file
 *  \brief External field ferm state and a creator
 */

#include "actions/ferm/fermstates/extfield_fermstate_w.h"
#include "actions/ferm/fermstates/extfield_aggregate_w.h"
#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

namespace Chroma
{

  /*! \ingroup fermstates */
  namespace CreateExtFieldFermStateEnv 
  { 
    CreateFermState<LatticeFermion,
		    multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> >* createFerm(XMLReader& xml, 
							      const std::string& path) 
    {
      
      return new CreateExtFieldFermState< LatticeFermion,
	multi1d<LatticeColorMatrix>, 
	multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml, path),
				      ExternalFieldEnv::reader(xml,path));
    }

    const std::string name = "EXTERNAL_FIELD_FERM_STATE";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories 
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheCreateFermStateFactory::Instance().registerObject(name, createFerm);

	success &= ExternalFieldEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
