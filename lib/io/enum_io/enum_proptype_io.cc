// -*- C++ -*-
// $Id: enum_proptype_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief PropType enum
 */
#include "enum_proptype_io.h"

#include <string>

namespace Chroma { 

  namespace PropTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = thePropTypeMap::Instance().registerPair(string("SCIDAC"), PROP_TYPE_SCIDAC );
      success &= thePropTypeMap::Instance().registerPair(string("SZIN"), PROP_TYPE_SZIN);
      success &= thePropTypeMap::Instance().registerPair( string("KYU"), PROP_TYPE_KYU );
      
      return success;
    }
    const string typeIDString = "PropType";
    bool registered = registerAll();
  };
  using namespace PropTypeEnv;

  //! Read a propagator type enum
  void read(XMLReader& xml_in,  const string& path, PropType& t) {
    thePropTypeMap::Instance().read(typeIDString,xml_in, path,t);
  }
  
  //! Write a propagator type enum
  void write(XMLWriter& xml_out, const string& path, const PropType& t) {
    thePropTypeMap::Instance().write(typeIDString,xml_out, path, t);
  }
};
