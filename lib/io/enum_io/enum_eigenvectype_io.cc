#include "enum_eigenvectype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace EigenVecTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theEigenVecTypeMap::Instance().registerPair(string("SCIDAC"), EVEC_TYPE_SCIDAC );
      success &=theEigenVecTypeMap::Instance().registerPair(string("SZIN"), EVEC_TYPE_SZIN);
      
      return success;
    }

    const bool registered = registerAll();
  };

  //! Read an eigenvectype enum
  void read(XMLReader& xml_in,  const string& path, EigenVecType& t) {
    theEigenVecTypeMap::Instance().read(xml_in, path,t);
  }
  
  //! Write an eigenvectype enum
  void write(XMLWriter& xml_out, const string& path, const EigenVecType& t) {
    theEigenVecTypeMap::Instance().write(xml_out, path, t);
  }
};
