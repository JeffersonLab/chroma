#include "chromabase.h"
#include "gaugebc.h"

#include "actions/gauge/gaugebc_simple.h"
#include "actions/gauge/gaugebc_factory.h"


using namespace QDP;
using namespace Chroma;
using namespace std;


namespace Chroma {

  namespace SimpleGaugeBCEnv { 

    //! Calllback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      QDPIO::cout << "Factory Callback: Creating SimpleGaugeBC " << endl;
      return new SimpleGaugeBC(SimpleGaugeBCParams(xml, path));
    }

    const std::string name = "GAUGEBC_SIMPLE";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);

  };

  SimpleGaugeBCParams::SimpleGaugeBCParams(XMLReader& xml, 
					   const std::string& path) {
      
    XMLReader paramtop(xml, path);
    
    try { 
      read(paramtop, "./boundary", boundary);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " << e << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Creating SimpleGaugeBCParams with boundary: ";
    QDPIO::cout << "Boundary.size = " << boundary.size() << endl;
    for(int i = 0; i < boundary.size();i++) { 
      QDPIO::cout << boundary[i] << endl;
    }
    QDPIO::cout << endl;
  }




  void read(XMLReader& xml, const std::string& path, SimpleGaugeBCParams& p) {
    SimpleGaugeBCParams tmp(xml, path);
    p=tmp;
  }
  
  void write(XMLWriter& xml, const std::string& path, const SimpleGaugeBCParams& p) { 
    push(xml, path);
    write(xml, "boundary", p.boundary);
    pop(xml);
  }

}; // End namespace Chroma 
