// $Id: gaugebc_aggregate.cc,v 3.3 2006-09-21 18:43:26 edwards Exp $
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
#include "actions/gauge/gaugebcs/schr_sf_zero_gaugebc.h"

namespace Chroma
{

  //! Name and registration
  namespace GaugeTypeGaugeBCEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= SimpleGaugeBCEnv::registerAll();
	success &= PeriodicGaugeBCEnv::registerAll();
	success &= SchrTrivialGaugeBCEnv::registerAll();
	success &= SchrNonPertGaugeBCEnv::registerAll();
	success &= SchrCouplingGaugeBCEnv::registerAll();
	success &= SchrChromoMagGaugeBCEnv::registerAll();
	success &= SchrDirichletGaugeBCEnv::registerAll();
	success &= SchrSFZeroGaugeBCEnv::registerAll();

	registered = true;
      }
      return success;
    }


    // Helper function for the GaugeAction readers
    Handle<GaugeBC< multi1d<LatticeColorMatrix>, 
		    multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
							   const std::string& path)
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

      Handle<GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	gbc(TheGaugeBCFactory::Instance().createObject(gaugebc,
						       top,
						       gaugebc_path));

      return gbc;
    }
  }

}
