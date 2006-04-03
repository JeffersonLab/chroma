// -*- C++ -*-
// $Id: enum_fermtype_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief FermType enum
 */

#include "enum_fermtype_io.h"

#include <string>

namespace Chroma { 

  namespace FermTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theFermTypeMap::Instance().registerPair(string("WILSON"), FERM_TYPE_WILSON );
      success &=theFermTypeMap::Instance().registerPair(string("STAGGERED"), FERM_TYPE_STAGGERED);
      
      return success;
    }
    const string typeIDString = "FermType";

    bool registered = registerAll();
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
