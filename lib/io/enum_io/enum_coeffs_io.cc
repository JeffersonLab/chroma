#include "enum_coeffs_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace CoeffTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theCoeffTypeMap::Instance().registerPair(string("ZOLOTAREV"), COEFF_TYPE_ZOLOTAREV );
      success &= theCoeffTypeMap::Instance().registerPair(string("TANH"), COEFF_TYPE_TANH);
      
      return success;
    }

    const bool registered = registerAll();
  };

  void read(XMLReader& xml_in,  const string& path, CoeffType& t) {
    theCoeffTypeMap::Instance().read(xml_in, path,t);
  }
  
  void write(XMLWriter& xml_out, const string& path, const CoeffType& t) {
    theCoeffTypeMap::Instance().write(xml_out, path, t);
  }
};
