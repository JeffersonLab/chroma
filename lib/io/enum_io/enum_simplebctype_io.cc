#include "enum_simplebctype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace SimpleBCTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theSimpleBCTypeMap::Instance().registerPair(string("ANTIPERIODIC"), BC_TYPE_ANTIPERIODIC );
      success &=theSimpleBCTypeMap::Instance().registerPair(string("DIRICHLET"), BC_TYPE_DIRICHLET);
      success &=theSimpleBCTypeMap::Instance().registerPair(string("PERIODIC"), BC_TYPE_PERIODIC);
      
      return success;
    }

    const bool registered = registerAll();
  };

  //! Read an simpleBC type enum
  void read(XMLReader& xml_in,  const string& path, SimpleBCType& t) {
    theSimpleBCTypeMap::Instance().read(xml_in, path,t);
  }
  
  //! Write an simpleBC type enum
  void write(XMLWriter& xml_out, const string& path, const SimpleBCType& t) {
    theSimpleBCTypeMap::Instance().write(xml_out, path, t);
  }
};
