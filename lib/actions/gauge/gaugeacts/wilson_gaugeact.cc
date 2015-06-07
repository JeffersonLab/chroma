/*! \file
 *  \brief Wilson gauge action
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma
{
 
  namespace WilsonGaugeActEnv 
  { 
    GaugeAction< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
											    const std::string& path) 
    {
      return new WilsonGaugeAct(CreateGaugeStateEnv::reader(xml, path), 
				WilsonGaugeActParams(xml, path));
    }

    const std::string name = "WILSON_GAUGEACT";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeActFactory::Instance().registerObject(name, createGaugeAct);
	registered = true;
      }
      return success;
    }
  }


  // Private initializer
  void
  WilsonGaugeAct::init(Handle< CreateGaugeState<P,Q> > cgs)
  {
    START_CODE();

    // Fold in normalizations and create action
    Real coeff = param.beta;
    plaq = new PlaqGaugeAct(cgs,coeff,param.aniso);
    
    END_CODE();
  } 

}

