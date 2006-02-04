// -*- C++ -*-
// $Id: enum_invtype_io.cc,v 2.1 2006-02-04 16:48:24 edwards Exp $
/*! \file
 * \brief Inverter type enum
 */
#include "enum_invtype_io.h"

#include <string>

namespace Chroma { 

  namespace InvTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theInvTypeMap::Instance().registerPair(string("CG_INVERTER"), CG_INVERTER );
      success &=theInvTypeMap::Instance().registerPair(string("MR_INVERTER"), MR_INVERTER);
      success &=theInvTypeMap::Instance().registerPair(string("BICG_INVERTER"), BICG_INVERTER);
      success &=theInvTypeMap::Instance().registerPair(string("SUMR_INVERTER"), SUMR_INVERTER);
      success &=theInvTypeMap::Instance().registerPair(string("REL_CG_INVERTER"), REL_CG_INVERTER);
      success &=theInvTypeMap::Instance().registerPair(string("REL_SUMR_INVERTER"), REL_SUMR_INVERTER);
      success &=theInvTypeMap::Instance().registerPair(string("REL_GMRESR_SUMR_INVERTER"), REL_GMRESR_SUMR_INVERTER);
      success &=theInvTypeMap::Instance().registerPair(string("REL_GMRESR_CG_INVERTER"), REL_GMRESR_CG_INVERTER);
      success &=theInvTypeMap::Instance().registerPair(string("BICGSTAB_INVERTER"), BICGSTAB_INVERTER);
      
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "InvType";
  };

  using namespace InvTypeEnv ;

  //! Read an InvType enum
  void read(XMLReader& xml_in,  const string& path, InvType& t) {
    theInvTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an InvType enum
  void write(XMLWriter& xml_out, const string& path, const InvType& t) {
    theInvTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
