// $Id: gauge_createstate_aggregate.cc,v 1.3 2006-09-21 19:42:07 edwards Exp $
/*! \file
 *  \brief All gauge create-state method
 */

#include "chromabase.h"

#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"

#include "actions/gauge/gaugestates/periodic_gaugestate.h"
#include "actions/gauge/gaugestates/simple_gaugestate.h"
#include "actions/gauge/gaugestates/stout_gaugestate.h"

#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"

namespace Chroma
{

  //! Registration aggregator
  namespace CreateGaugeStateEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Register all gauge BCs
	success &= GaugeTypeGaugeBCEnv::registerAll();

	// Register all gauge states
	success &= CreatePeriodicGaugeStateEnv::registerAll();
	success &= CreateSimpleGaugeStateEnv::registerAll();
	success &= CreateStoutGaugeStateEnv::registerAll();

	registered = true;
      }
      return success;
    }


    // Helper function for the GaugeAction readers
    Handle< CreateGaugeState< multi1d<LatticeColorMatrix>, 
			      multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								     const std::string& path)
    {
      XMLReader top(xml_in, path);

      std::string gaugestate;
      std::string gaugestate_path;
      if (top.count("GaugeState") != 0)
      {
	gaugestate_path = "GaugeState";
	read(top, gaugestate_path + "/Name", gaugestate);
      }
      else
      {
//	QDPIO::cerr << "Error: GaugeState group not found" << endl;
//	QDP_abort(1);
	
	gaugestate_path = ".";
	gaugestate = Chroma::CreateSimpleGaugeStateEnv::name;
      }

      Handle< CreateGaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	cgs(TheCreateGaugeStateFactory::Instance().createObject(gaugestate,
								top,
								gaugestate_path));

      return cgs;
    }


    // Returns a no-smearing group
    GroupXML_t   nullXMLGroup()
    {
      GroupXML_t nope;
      nope.id = CreatePeriodicGaugeStateEnv::name;
      nope.path = "GaugeState";

      XMLBufferWriter xml_tmp;
      push(xml_tmp, "GaugeState");
      write(xml_tmp, "Name", nope.id);
      pop(xml_tmp);

      nope.xml = xml_tmp.str();

      return nope;
    }

  }
}
