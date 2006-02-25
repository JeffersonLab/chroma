/*! \file
 *  \brief Periodic gauge boundary conditions
 */

#include "chromabase.h"
#include "gaugebc.h"

#include "actions/gauge/gaugebcs/gaugebc_periodic.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace PeriodicGaugeBCEnv 
  { 

    //! Callback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      QDPIO::cout << "GaugeBC Callback: Creating PeriodicGaugeBC " << endl;
      return new PeriodicGaugeBC();
    }

    const std::string name = "PERIODIC_GAUGEBC";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);

  };

}; // End namespace Chroma 
