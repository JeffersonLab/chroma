// $Id: stout_gaugestate.cc,v 3.1 2006-09-07 01:26:43 bjoo Exp $
/*! \file
 *  \brief Simple gauge state and a creator
 */

#include "actions/gauge/gaugeacts/stout_gaugestate.h"
#include "actions/gauge/gaugeacts/gauge_createstate_factory.h"
#include "actions/gauge/gaugeacts/gauge_createstate_aggregate.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"

namespace Chroma
{
 
  namespace CreateStoutGaugeStateEnv 
  { 
    CreateGaugeState<multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> >* createCreator(XMLReader& xml, 
								  const std::string& path) 
    {
      return new CreateStoutGaugeState<multi1d<LatticeColorMatrix>, 
	multi1d<LatticeColorMatrix> >(GaugeTypeGaugeBCEnv::reader(xml, path),
				      StoutFermStateParams(xml, path));
    }

    const std::string name = "STOUT_GAUGE_STATE";
    const bool registered = TheCreateGaugeStateFactory::Instance().registerObject(name, 
										  createCreator);
  }

}

