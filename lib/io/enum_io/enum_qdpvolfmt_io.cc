#include "enum_qdpvolfmt_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace QDPVolfmtEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theQDPVolfmtMap::Instance().registerPair(string("SINGLEFILE"), QDPIO_SINGLEFILE );
      success &=theQDPVolfmtMap::Instance().registerPair(string("MULTIFILE"), QDPIO_MULTIFILE);
      
      return success;
    }

    const bool registered = registerAll();
  };

  //! Read a QDP volume format type
  void read(XMLReader& xml_in,  const string& path, QDP_volfmt_t& t) {
    theQDPVolfmtMap::Instance().read(xml_in, path,t);
  }
  
  //! Write a QDP volume format type
  void write(XMLWriter& xml_out, const string& path, const QDP_volfmt_t& t) {
    theQDPVolfmtMap::Instance().write(xml_out, path, t);
  }
};
