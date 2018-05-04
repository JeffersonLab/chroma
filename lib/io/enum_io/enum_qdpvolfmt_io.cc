// -*- C++ -*-
/*! \file
 * \brief QDP volfmt enum
 */
#include "enum_qdpvolfmt_io.h"

#include <string>

namespace Chroma { 

  namespace QDPVolfmtEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theQDPVolfmtMap::Instance().registerPair(std::string("SINGLEFILE"), QDPIO_SINGLEFILE );
      success &=theQDPVolfmtMap::Instance().registerPair(std::string("MULTIFILE"), QDPIO_MULTIFILE);
      success &=theQDPVolfmtMap::Instance().registerPair(std::string("PARTFILE"), QDPIO_PARTFILE);
      
      return success;
    }

    const std::string typeIDString = "QDP_volfmt_t";

    bool registered = registerAll();
  };
  using namespace QDPVolfmtEnv;

  //! Read a QDP volume format type
  void read(XMLReader& xml_in,  const std::string& path, QDP_volfmt_t& t) {
    theQDPVolfmtMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write a QDP volume format type
  void write(XMLWriter& xml_out, const std::string& path, const QDP_volfmt_t& t) {
    theQDPVolfmtMap::Instance().write(typeIDString,xml_out, path, t);
  }
};
