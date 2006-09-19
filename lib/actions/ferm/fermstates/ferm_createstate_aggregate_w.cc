// $Id: ferm_createstate_aggregate_w.cc,v 1.1 2006-09-19 17:53:36 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#include "chromabase.h"

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"

#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermstates/periodic_fermstate_w.h"
#include "actions/ferm/fermstates/simple_fermstate_w.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace CreateFermStateEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // All ferm bcs
      success &= WilsonTypeFermBCEnv::registered;

      // All fermstates
      success &= CreatePeriodicFermStateEnv::registered;
      success &= CreateSimpleFermStateEnv::registered;
      success &= CreateStoutFermStateEnv::registered;
      success &= CreateSLICFermStateEnv::registered;

      return success;
    }

    const bool registered = registerAll();


    // Returns a periodic group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;
      nope.id = CreatePeriodicFermStateEnv::name;

      XMLBufferWriter xml_tmp;
      push(xml_tmp, "FermState");
      write(xml_tmp, "Name", nope.id);
      pop(xml_tmp);

      nope.xml = xml_tmp.str();

      return nope;
    }

  }
}
