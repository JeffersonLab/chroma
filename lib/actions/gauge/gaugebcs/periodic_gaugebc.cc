/*! \file
 *  \brief Periodic gauge boundary conditions
 */

#include "chromabase.h"
#include "gaugebc.h"

#include "actions/gauge/gaugebcs/gaugebc_factory.h"
#include "actions/gauge/gaugebcs/periodic_gaugebc.h"

namespace Chroma 
{

  namespace PeriodicGaugeBCEnv 
  { 

    //! Callback function to register with the factory
    GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeBC(XMLReader& xml, 
										       const string& path)
    {
      QDPIO::cout << "GaugeBC Callback: Creating PeriodicGaugeBC " << endl;
      return new PeriodicGaugeBC();
    }

    const std::string name = "PERIODIC_GAUGEBC";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeBCFactory::Instance().registerObject(name, createGaugeBC);
	registered = true;
      }
      return success;
    }
  }

} // End namespace Chroma 
