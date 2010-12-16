// $Id: ferm_createstate_aggregate_w.cc,v 1.5 2007-11-01 20:56:51 kostas Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#include "chromabase.h"

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"

#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermstates/periodic_fermstate_w.h"
#include "actions/ferm/fermstates/simple_fermstate_w.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "actions/ferm/fermstates/hex_fermstate_w.h"
#include "actions/ferm/fermstates/extfield_fermstate_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace CreateFermStateEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// All ferm bcs
	success &= WilsonTypeFermBCEnv::registerAll();

	// All fermstates
	success &= CreatePeriodicFermStateEnv::registerAll();
	success &= CreateSimpleFermStateEnv::registerAll();
	success &= CreateStoutFermStateEnv::registerAll();
	success &= CreateHexFermStateEnv::registerAll();
	success &= CreateSLICFermStateEnv::registerAll();
	success &= CreateExtFieldFermStateEnv::registerAll();

	registered = true;
      }
      return success;
    }


    // Returns a periodic group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;
      nope.id = CreatePeriodicFermStateEnv::name;
      nope.path = "FermState";

      XMLBufferWriter xml_tmp;
      push(xml_tmp, "FermState");
      write(xml_tmp, "Name", nope.id);
      pop(xml_tmp);

      nope.xml = xml_tmp.str();

      return nope;
    }

  }
}
