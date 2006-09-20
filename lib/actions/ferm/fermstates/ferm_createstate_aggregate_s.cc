// $Id: ferm_createstate_aggregate_s.cc,v 1.2 2006-09-20 20:31:41 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#include "chromabase.h"

#include "actions/ferm/fermstates/ferm_createstate_aggregate_s.h"
#include "actions/ferm/fermstates/periodic_fermstate_s.h"
#include "actions/ferm/fermstates/simple_fermstate_s.h"

namespace Chroma
{

  //! Registration aggregator
  namespace StaggeredCreateFermStateEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= StaggeredCreatePeriodicFermStateEnv::registerAll();
	success &= StaggeredCreateSimpleFermStateEnv::registerAll();
	registered = true;
      }
      return success;
    }


    // Returns a periodic group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;
      nope.id = StaggeredCreatePeriodicFermStateEnv::name;

      XMLBufferWriter xml_tmp;
      push(xml_tmp, "FermState");
      write(xml_tmp, "Name", nope.id);
      pop(xml_tmp);

      nope.xml = xml_tmp.str();

      return nope;
    }

  }

}
