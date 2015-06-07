// -*- C++ -*-
/*! \file
 * \brief Simple BC type enum
 */
#include "enum_simplebctype_io.h"

#include <string>

namespace Chroma { 

  namespace SimpleBCTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theSimpleBCTypeMap::Instance().registerPair(std::string("ANTIPERIODIC"), BC_TYPE_ANTIPERIODIC );
      success &=theSimpleBCTypeMap::Instance().registerPair(std::string("DIRICHLET"), BC_TYPE_DIRICHLET);
      success &=theSimpleBCTypeMap::Instance().registerPair(std::string("PERIODIC"), BC_TYPE_PERIODIC);
      
      return success;
    }

    // Boilerplate stuff from here on
    const std::string typeIDString = "SimpleBCType";
    bool registered = registerAll();
  };
  using namespace SimpleBCTypeEnv;

  //! Read an simpleBC type enum
  void read(XMLReader& xml_in,  const std::string& path, SimpleBCType& t) {
    theSimpleBCTypeMap::Instance().read(typeIDString,xml_in, path,t);
  }
  
  //! Write an simpleBC type enum
  void write(XMLWriter& xml_out, const std::string& path, const SimpleBCType& t) {
    theSimpleBCTypeMap::Instance().write(typeIDString,xml_out, path, t);
  }
};
