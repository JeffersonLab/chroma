// -*- C++ -*-
// $Id: enum_plusminus_io.cc,v 1.1 2006-05-23 18:11:34 edwards Exp $
/*! \file
 * \brief Inverter type enum
 */
#include "enum_plusminus_io.h"

#include <string>

namespace Chroma 
{ 

  namespace PlusMinusEnv 
  { 

    bool registerAll(void) 
    {
      bool success = true; 
      success &= thePlusMinusMap::Instance().registerPair(std::string("PLUS"), PLUS);
      success &= thePlusMinusMap::Instance().registerPair(std::string("MINUS"), MINUS);
      return success;
    }

    bool registered = registerAll();
    const std::string typeIDString = "PlusMinus";
  };

  using namespace PlusMinusEnv ;

  //! Read an PlusMinus enum
  void read(XMLReader& xml_in,  const std::string& path, PlusMinus& t) {
    thePlusMinusMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an PlusMinus enum
  void write(XMLWriter& xml_out, const std::string& path, const PlusMinus& t) {
    thePlusMinusMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
