#include "enum_stochsrc_io.h"

#include <string>

namespace Chroma { 

  namespace StochSrcEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theStochSrc::Instance().registerPair(string("Z2NOISE"), Z2NOISE );
      success &=theStochSrc::Instance().registerPair(string("GAUSSIAN"),GAUSSIAN );
      
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "StochSrc";
  };

  using namespace StochSrcEnv ;

  //! Read an StochSrc enum
  void read(XMLReader& xml_in,  const string& path, VolSrc& t) {
    theStochSrc::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an StochSrc enum
  void write(XMLWriter& xml_out, const string& path, const VolSrc& t) {
    theStochSrc::Instance().write(typeIDString, xml_out, path, t);
  }
};
