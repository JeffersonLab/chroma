// -*- C++ -*-
/*! \file
 * \brief Heatbath enum
 */

#include "enum_heatbathtype_io.h"

#include <string>

namespace Chroma { 

  namespace HeatbathTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theHeatbathTypeMap::Instance().registerPair(std::string("KPHB"), HEATBATH_TYPE_KPHB );
      success &=theHeatbathTypeMap::Instance().registerPair(std::string("CrHB"), HEATBATH_TYPE_CrHB);
      
      return success;
    }

    const std::string typeIDString ="HeatbathType";
    bool registered = registerAll();
  };

  using namespace HeatbathTypeEnv;
  //! Read an HeatbathType enum
  void read(XMLReader& xml_in,  const std::string& path, HeatbathType& t) {
    theHeatbathTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an HeatbathType enum
  void write(XMLWriter& xml_out, const std::string& path, const HeatbathType& t) {
    theHeatbathTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
