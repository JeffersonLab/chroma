#include "enum_fermtype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace FermTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theFermTypeMap::Instance().registerPair(string("WILSON"), FERM_TYPE_WILSON );
      success &=theFermTypeMap::Instance().registerPair(string("STAGGERED"), FERM_TYPE_STAGGERED);
      
      return success;
    }
    const string typeIDString = "FermType";

    const bool registered = registerAll();
  };

  using namespace FermTypeEnv;

  //! Read an fermion type enum
  void read(XMLReader& xml_in,  const string& path, FermType& t) {
    theFermTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an fermion type enum
  void write(XMLWriter& xml_out, const string& path, const FermType& t) {
    theFermTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
