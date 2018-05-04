// -*- C++ -*-
/*! \file
 * \brief Eigenstd::vector type enum
 */

#include "enum_eigenvectype_io.h"

#include <string>

namespace Chroma { 

  namespace EigenVecTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theEigenVecTypeMap::Instance().registerPair(std::string("SCIDAC"), EVEC_TYPE_SCIDAC );
      success &=theEigenVecTypeMap::Instance().registerPair(std::string("SZIN"), EVEC_TYPE_SZIN);
      
      return success;
    }
    const std::string typeIDString = "EigenVecType" ;
    bool registered = registerAll();
  };
  using namespace EigenVecTypeEnv;

  //! Read an eigenvectype enum
  void read(XMLReader& xml_in,  const std::string& path, EigenVecType& t) {
    theEigenVecTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an eigenvectype enum
  void write(XMLWriter& xml_out, const std::string& path, const EigenVecType& t) {
    theEigenVecTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
