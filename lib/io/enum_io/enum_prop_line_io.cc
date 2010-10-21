// -*- C++ -*-
/*! \file
 * \brief Type of lines for distillution
 */

#include "enum_prop_line_io.h"

#include <string>

namespace Chroma { 

  namespace PropLineTypeEnv { 

    bool registerAll(void) 
    {
      bool success = true; 
      success &= thePropLineTypeMap::Instance().registerPair(string("CONN"), PROP_LINE_TYPE_CONN);
      success &= thePropLineTypeMap::Instance().registerPair(string("DISC"), PROP_LINE_TYPE_DISC);
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "PropLineType";
  };
  using namespace PropLineTypeEnv;

  //! Reader
  void read(XMLReader& xml_in, const string& path, PropLineType& t) {
    thePropLineTypeMap::Instance().read(typeIDString, xml_in, path, t);
  }
  
  //! Writer
  void write(XMLWriter& xml_out, const string& path, const PropLineType& t) {
    thePropLineTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }

  //! Reader
  void read(BinaryReader& bin_in, PropLineType& t) {
    thePropLineTypeMap::Instance().read(typeIDString, bin_in, t);
  }
  
  //! Writer
  void write(BinaryWriter& bin_out, const PropLineType& t) {
    thePropLineTypeMap::Instance().write(typeIDString, bin_out, t);
  }
};
