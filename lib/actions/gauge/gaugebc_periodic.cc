#include "chromabase.h"
#include "gaugebc.h"

#include "actions/gauge/gaugebc_periodic.h"
#include "actions/gauge/gaugebc_factory.h"


using namespace QDP;
using namespace Chroma;
using namespace std;


namespace Chroma {

  namespace PeriodicGaugeBCEnv { 

    //! Callback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      QDPIO::cout << "GaugeBC Callback: Creating PeriodicGaugeBC " << endl;
      return new PeriodicGaugeBC();
    }

    const std::string name = "GAUGEBC_PERIODIC";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);

  };

}; // End namespace Chroma 
