// -*- C++ -*-
// $Id: enum_stochsrc_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief Stochastic source enum
 */
#include "enum_stochsrc_io.h"

#include <string>

namespace Chroma { 

  namespace StochSrcEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theStochSrc::Instance().registerPair(string("Z2NOISE"), Z2NOISE );
      success &=theStochSrc::Instance().registerPair(string("GAUSSIAN"),GAUSSIAN );


      success &=theStochSrc::Instance().registerPair(string("T_DILUTE_GAUSS"),T_DILUTE_GAUSS );
      success &=theStochSrc::Instance().registerPair(string("C_DILUTE_GAUSS"),C_DILUTE_GAUSS );
      success &=theStochSrc::Instance().registerPair(string("P_DILUTE_GAUSS"),P_DILUTE_GAUSS );
      success &=theStochSrc::Instance().registerPair(string("CT_DILUTE_GAUSS"),CT_DILUTE_GAUSS );
      success &=theStochSrc::Instance().registerPair(string("CP_DILUTE_GAUSS"),CP_DILUTE_GAUSS );
      success &=theStochSrc::Instance().registerPair(string("PT_DILUTE_GAUSS"),PT_DILUTE_GAUSS );
      success &=theStochSrc::Instance().registerPair(string("MOD_T_DILUTE_GAUSS"),  MOD_T_DILUTE_GAUSS);
      success &=theStochSrc::Instance().registerPair(string("CORNER_DILUTE_GAUSS"),CORNER_DILUTE_GAUSS );
      success &=theStochSrc::Instance().registerPair(string("COR_DBL_T_DILUTE_GAUSS"), COR_DBL_T_DILUTE_GAUSS);
      success &=theStochSrc::Instance().registerPair(string("COR_MOD_DBL_T_DILUTE_GAUSS"), COR_MOD_DBL_T_DILUTE_GAUSS);
      success &=theStochSrc::Instance().registerPair(string("C_MOD_T_DILUTE_GAUSS"),C_MOD_T_DILUTE_GAUSS );

      
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "StochSrc";
  };

  using namespace StochSrcEnv ;

  //! Read an StochSrc enum
  void read(XMLReader& xml_in,  const string& path, VolSrc& t) {
    theStochSrc::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an StochSrc enum
  void write(XMLWriter& xml_out, const string& path, const VolSrc& t) {
    theStochSrc::Instance().write(typeIDString, xml_out, path, t);
  }
};
