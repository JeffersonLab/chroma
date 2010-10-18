// -*- C++ -*-
/*! \file
 * \brief Type of lines for distillution
 */

#include "enum_prop_dist_io.h"

#include <string>

namespace Chroma { 

  namespace PropDistTypeEnv { 

    bool registerAll(void) 
    {
      bool success = true; 
      success &= thePropDistTypeMap::Instance().registerPair(string("SINGLE-SOURCE"), PROP_DIST_TYPE_SINGLE_SOURCE);
      success &= thePropDistTypeMap::Instance().registerPair(string("MULTI-SOURCE"), PROP_DIST_TYPE_MULTI_SOURCE);
      success &= thePropDistTypeMap::Instance().registerPair(string("SINGLE-SOLUTION"), PROP_DIST_TYPE_SINGLE_SOLUTION);
      success &= thePropDistTypeMap::Instance().registerPair(string("MULTI-SOLUTION"), PROP_DIST_TYPE_MULTI_SOLUTION);
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "PropDistType";
  };
  using namespace PropDistTypeEnv;

  //! Reader
  void read(XMLReader& xml_in, const string& path, PropDistType& t) {
    thePropDistTypeMap::Instance().read(typeIDString, xml_in, path, t);
  }
  
  //! Writer
  void write(XMLWriter& xml_out, const string& path, const PropDistType& t) {
    thePropDistTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }

  //! Reader
  void read(BinaryReader& bin_in, PropDistType& t) {
    thePropDistTypeMap::Instance().read(typeIDString, bin_in, t);
  }
  
  //! Writer
  void write(BinaryWriter& bin_out, const PropDistType& t) {
    thePropDistTypeMap::Instance().write(typeIDString, bin_out, t);
  }
};
