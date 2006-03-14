// $Id: gaugebc_aggregate.cc,v 2.3 2006-03-14 04:49:54 edwards Exp $
/*! \file
 *  \brief Gauge boundary condition aggregator
 */

#include "actions/gauge/gaugebcs/gaugebc_factory.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"
#include "actions/gauge/gaugebcs/simple_gaugebc.h"
#include "actions/gauge/gaugebcs/periodic_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_triv_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_nonpert_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_coupling_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_chromomag_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_dirich_gaugebc.h"

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
      success &= SchrTrivialGaugeBCEnv::registered;
      success &= SchrNonPertGaugeBCEnv::registered;
      success &= SchrCouplingGaugeBCEnv::registered;
      success &= SchrChromoMagGaugeBCEnv::registered;
      success &= SchrDirichletGaugeBCEnv::registered;
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
