// -*- C++ -*-
/*! \file
 * \brief Coeffs enum
 */

#include "enum_coeffs_io.h"

#include <string>

namespace Chroma { 

  namespace CoeffTypeEnv { 

    bool registerAll(void) 
    {
      bool success = true; 
      success &= theCoeffTypeMap::Instance().registerPair(std::string("ZOLOTAREV"), COEFF_TYPE_ZOLOTAREV );
      success &= theCoeffTypeMap::Instance().registerPair(std::string("TANH"), COEFF_TYPE_TANH);
      success &= theCoeffTypeMap::Instance().registerPair(std::string("TANH_UNSCALED"), COEFF_TYPE_TANH_UNSCALED);
      return success;
    }

    bool registered = registerAll();
    const std::string typeIDString = "CoeffType";
  };
  using namespace CoeffTypeEnv;

  //! read an approximation coefficient type enum
  void read(XMLReader& xml_in,  const std::string& path, CoeffType& t) {
    theCoeffTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! write an approximation coefficient type enum
  void write(XMLWriter& xml_out, const std::string& path, const CoeffType& t) {
    theCoeffTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
