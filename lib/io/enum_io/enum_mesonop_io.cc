// -*- C++ -*-
// $Id: enum_mesonop_io.cc,v 1.1 2008-06-04 03:16:48 edwards Exp $
/*! \file
 * \brief Type of contractions for stochastic meson operators
 */

#include "enum_mesonop_io.h"

#include <string>

namespace Chroma { 

  namespace MesonOpTypeEnv { 

    bool registerAll(void) 
    {
      bool success = true; 
      success &= theMesonOpTypeMap::Instance().registerPair(string("SOURCE-SOURCE"), MESON_OP_TYPE_SOURCE_SOURCE);
      success &= theMesonOpTypeMap::Instance().registerPair(string("SOURCE-SOLUTION"), MESON_OP_TYPE_SOURCE_SOLUTION);
      success &= theMesonOpTypeMap::Instance().registerPair(string("SOLUTION-SOURCE "), MESON_OP_TYPE_SOLUTION_SOURCE);
      success &= theMesonOpTypeMap::Instance().registerPair(string("SOLUTION-SOLUTION"), MESON_OP_TYPE_SOLUTION_SOLUTION);
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "MesonOpType";
  };
  using namespace MesonOpTypeEnv;

  //! read an approximation coefficient type enum
  void read(XMLReader& xml_in,  const string& path, MesonOpType& t) {
    theMesonOpTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! write an approximation coefficient type enum
  void write(XMLWriter& xml_out, const string& path, const MesonOpType& t) {
    theMesonOpTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
