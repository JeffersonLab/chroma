#include "enum_heatbathtype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace HeatbathTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theHeatbathTypeMap::Instance().registerPair(string("KPHB"), HEATBATH_TYPE_KPHB );
      success &=theHeatbathTypeMap::Instance().registerPair(string("CrHB"), HEATBATH_TYPE_CrHB);
      
      return success;
    }

    const bool registered = registerAll();
  };

  //! Read an fermion type enum
  void read(XMLReader& xml_in,  const string& path, HeatbathType& t) {
    theHeatbathTypeMap::Instance().read(xml_in, path,t);
  }
  
  //! Write an fermion type enum
  void write(XMLWriter& xml_out, const string& path, const HeatbathType& t) {
    theHeatbathTypeMap::Instance().write(xml_out, path, t);
  }
};
