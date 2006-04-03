// -*- C++ -*-
// $Id: enum_coeffs_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
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
      success &= theCoeffTypeMap::Instance().registerPair(string("ZOLOTAREV"), COEFF_TYPE_ZOLOTAREV );
      success &= theCoeffTypeMap::Instance().registerPair(string("TANH"), COEFF_TYPE_TANH);
      success &= theCoeffTypeMap::Instance().registerPair(string("TANH_UNSCALED"), COEFF_TYPE_TANH_UNSCALED);
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "CoeffType";
  };
  using namespace CoeffTypeEnv;

  //! read an approximation coefficient type enum
  void read(XMLReader& xml_in,  const string& path, CoeffType& t) {
    theCoeffTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! write an approximation coefficient type enum
  void write(XMLWriter& xml_out, const string& path, const CoeffType& t) {
    theCoeffTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
