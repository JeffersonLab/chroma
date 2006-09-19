// $Id: extfield_fermstate_w.cc,v 1.1 2006-09-19 17:53:36 edwards Exp $
/*! \file
 *  \brief External field ferm state and a creator
 */

#if 0

#error "TURNED OFF FOR THE MOMENT - NEED TO FIX CREATEEXTFIELDFermState CONSTRUCTOR CALL"


#include "actions/ferm/fermstates/extfield_fermstate_w.h"
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
	multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml, 
								  path));
    }

    const std::string name = "EXTERNAL_FIELD_FERM_STATE";

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


#endif
