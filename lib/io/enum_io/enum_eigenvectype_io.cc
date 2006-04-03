// -*- C++ -*-
// $Id: enum_eigenvectype_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief Eigenvector type enum
 */

#include "enum_eigenvectype_io.h"

#include <string>

namespace Chroma { 

  namespace EigenVecTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theEigenVecTypeMap::Instance().registerPair(string("SCIDAC"), EVEC_TYPE_SCIDAC );
      success &=theEigenVecTypeMap::Instance().registerPair(string("SZIN"), EVEC_TYPE_SZIN);
      
      return success;
    }
    const string typeIDString = "EigenVecType" ;
    bool registered = registerAll();
  };
  using namespace EigenVecTypeEnv;

  //! Read an eigenvectype enum
  void read(XMLReader& xml_in,  const string& path, EigenVecType& t) {
    theEigenVecTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an eigenvectype enum
  void write(XMLWriter& xml_out, const string& path, const EigenVecType& t) {
    theEigenVecTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
