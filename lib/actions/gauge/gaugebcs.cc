// $Id: gaugebcs.cc,v 1.1 2005-01-12 04:44:53 edwards Exp $
/*! \file
 *  \brief All gauge BC
 */

#include "actions/gauge/gaugebcs.h"

#warning "OTHER STUFF"
#include "actions/gauge/gaugebc_factory.h"
#include "actions/gauge/gaugebc_simple.h"
#include "actions/gauge/gaugebc_periodic.h"

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
