#include "enum_gaugeacttype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace GaugeActTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theGaugeActTypeMap::Instance().registerPair(string("WILSON"), GAUGE_ACT_TYPE_WILSON );
      
      return success;
    }

    const bool registered = registerAll();
  };

  //! Read an GaugeActType enum
  void read(XMLReader& xml_in,  const string& path, GaugeActType& t) {
    theGaugeActTypeMap::Instance().read(xml_in, path,t);
  }
  
  //! Write an GaugeActType enum
  void write(XMLWriter& xml_out, const string& path, const GaugeActType& t) {
    theGaugeActTypeMap::Instance().write(xml_out, path, t);
  }
};
