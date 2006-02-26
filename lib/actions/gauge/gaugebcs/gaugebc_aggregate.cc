// $Id: gaugebc_aggregate.cc,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Gauge boundary condition aggregator
 */

#include "actions/gauge/gaugebcs/gaugebc_factory.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"
#include "actions/gauge/gaugebcs/simple_gaugebc.h"
#include "actions/gauge/gaugebcs/periodic_gaugebc.h"
#include "actions/gauge/gaugebcs/schr1link_gaugebc.h"
#include "actions/gauge/gaugebcs/schr2link_gaugebc.h"

namespace Chroma
{

  //! Name and registration
  namespace GaugeTypeGaugeBCEnv
  {
    bool registerAll(void) 
    {
      bool success = true; 
      success &= SimpleGaugeBCEnv::registered;
      success &= PeriodicGaugeBCEnv::registered;
      success &= Schr1LinkGaugeBCEnv::registered;
      success &= Schr2LinkGaugeBCEnv::registered;
      return success;
    }

    const bool registered = registerAll();

    // Helper function for the GaugeAction readers
    Handle<GaugeBC> reader(XMLReader& xml_in, const std::string& path)
    {
      XMLReader top(xml_in, path);

      bool success = registered;  // make sure all codes loaded

      std::string gaugebc;
      std::string gaugebc_path;
      if (top.count("GaugeBC") != 0)
      {
	gaugebc_path = "GaugeBC";
	read(top, gaugebc_path + "/Name", gaugebc);
      }
      else
      {
	QDPIO::cerr << "Error: GaugeBC group not found" << endl;
	QDP_abort(1);
      }

      Handle<GaugeBC> 
	gbc(TheGaugeBCFactory::Instance().createObject(gaugebc,
						       top,
						       gaugebc_path));

      return gbc;
    }
  }

}
