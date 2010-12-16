// -*- C++ -*-
// $Id: enum_cfgtype_io.cc,v 3.1 2006-04-16 03:06:49 edwards Exp $
/*! \file
 * \brief CfgType enum
 */

#include "enum_cfgtype_io.h"
#include <string>


namespace Chroma { 

  namespace CfgTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theCfgTypeMap::Instance().registerPair(string("MILC"), CFG_TYPE_MILC );
      success &= theCfgTypeMap::Instance().registerPair(string("NERSC"), CFG_TYPE_NERSC);
      success &= theCfgTypeMap::Instance().registerPair( string("SCIDAC"), CFG_TYPE_SCIDAC );
      success &= theCfgTypeMap::Instance().registerPair( string("SZIN" ), CFG_TYPE_SZIN );
      success &= theCfgTypeMap::Instance().registerPair( string("SZINQIO"), CFG_TYPE_SZINQIO );
      
      success &= theCfgTypeMap::Instance().registerPair( string("KYU"), CFG_TYPE_KYU );
      success &= theCfgTypeMap::Instance().registerPair( string("WUP"), CFG_TYPE_WUPP);


      success &= theCfgTypeMap::Instance().registerPair( string("DISORDERED"), CFG_TYPE_DISORDERED );
      success &= theCfgTypeMap::Instance().registerPair( string("UNIT"), CFG_TYPE_UNIT );
      success &= theCfgTypeMap::Instance().registerPair( string("CPPACS"), CFG_TYPE_CPPACS );
      success &= theCfgTypeMap::Instance().registerPair( string("WEAK_FIELD"), CFG_TYPE_WEAK_FIELD);
      success &= theCfgTypeMap::Instance().registerPair( string("CLASSICAL_SF"), CFG_TYPE_CLASSICAL_SF);
      
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "CfgType";
  };

  using namespace CfgTypeEnv;

  void read(XMLReader& xml_in, const string& path, CfgType& t) {
    theCfgTypeMap::Instance().read(typeIDString, xml_in,path,t);
  }
  
  void write(XMLWriter& xml_out, const string& path, const CfgType& t) {
    theCfgTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
