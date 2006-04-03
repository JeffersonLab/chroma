// -*- C++ -*-
// $Id: enum_qdpvolfmt_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
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
      success = theQDPVolfmtMap::Instance().registerPair(string("SINGLEFILE"), QDPIO_SINGLEFILE );
      success &=theQDPVolfmtMap::Instance().registerPair(string("MULTIFILE"), QDPIO_MULTIFILE);
      success &=theQDPVolfmtMap::Instance().registerPair(string("PARTFILE"), QDPIO_PARTFILE);
      
      return success;
    }

    const string typeIDString = "QDP_volfmt_t";

    bool registered = registerAll();
  };
  using namespace QDPVolfmtEnv;

  //! Read a QDP volume format type
  void read(XMLReader& xml_in,  const string& path, QDP_volfmt_t& t) {
    theQDPVolfmtMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write a QDP volume format type
  void write(XMLWriter& xml_out, const string& path, const QDP_volfmt_t& t) {
    theQDPVolfmtMap::Instance().write(typeIDString,xml_out, path, t);
  }
};
