// $Id: periodic_gaugestate.cc,v 1.3 2009-04-17 02:05:36 bjoo Exp $
/*! \file
 *  \brief Periodic gauge state and a creator
 */

#include "actions/gauge/gaugestates/periodic_gaugestate.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"

namespace Chroma
{

  /*! \ingroup gaugestates */
  namespace CreatePeriodicGaugeStateEnv 
  { 
    CreateGaugeState<multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> >* createCreator(XMLReader& xml, 
								  const std::string& path) 
    {
      return new CreatePeriodicGaugeState<multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> > ();
    }

    const std::string name = "PERIODIC_GAUGE_STATE";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheCreateGaugeStateFactory::Instance().registerObject(name, createCreator);
	registered = true;
      }
      return success;
    }
  }

}

