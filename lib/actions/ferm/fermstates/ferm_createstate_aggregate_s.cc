// $Id: ferm_createstate_aggregate_s.cc,v 1.1 2006-09-19 17:53:36 edwards Exp $
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
    bool registerAll() 
    {
      bool success = true;

      success &= StaggeredCreatePeriodicFermStateEnv::registered;
      success &= StaggeredCreateSimpleFermStateEnv::registered;

      return success;
    }

    const bool registered = registerAll();


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
