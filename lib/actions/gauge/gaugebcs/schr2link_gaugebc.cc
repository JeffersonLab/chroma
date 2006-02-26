// $Id: schr2link_gaugebc.cc,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Schroedinger functional 2-link gauge boundary conditions
 */

#include "actions/gauge/gaugebcs/schr2link_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace Schr2LinkGaugeBCEnv 
  { 
    //! Callback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      QDPIO::cout << "GaugeBC Callback: Creating SchrGaugeBC " << endl;
      return new Schr2LinkGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_2LINK_GAUGEBC";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);
  }

}
