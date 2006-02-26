// $Id: schr1link_gaugebc.cc,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Schroedinger functional 1-link gauge boundary conditions
 */
#include "chromabase.h"
#include "gaugebc.h"

#include "actions/gauge/gaugebcs/schr1link_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace Schr1LinkGaugeBCEnv 
  { 
    //! Callback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      QDPIO::cout << "GaugeBC Callback: Creating SchrGaugeBC " << endl;
      return new Schr1LinkGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_1LINK_GAUGEBC";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);
  }

}
