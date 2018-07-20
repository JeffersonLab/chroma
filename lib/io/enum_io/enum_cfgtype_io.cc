// -*- C++ -*-
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
      success = theCfgTypeMap::Instance().registerPair(std::string("MILC"), CFG_TYPE_MILC );
      success &= theCfgTypeMap::Instance().registerPair(std::string("NERSC"), CFG_TYPE_NERSC);
      success &= theCfgTypeMap::Instance().registerPair( std::string("SCIDAC"), CFG_TYPE_SCIDAC );
      success &= theCfgTypeMap::Instance().registerPair( std::string("SZIN" ), CFG_TYPE_SZIN );
      success &= theCfgTypeMap::Instance().registerPair( std::string("SZINQIO"), CFG_TYPE_SZINQIO );
      
      success &= theCfgTypeMap::Instance().registerPair( std::string("KYU"), CFG_TYPE_KYU );
      success &= theCfgTypeMap::Instance().registerPair( std::string("WUP"), CFG_TYPE_WUPP);
      success &= theCfgTypeMap::Instance().registerPair( std::string("DISORDERED"), CFG_TYPE_DISORDERED );
      success &= theCfgTypeMap::Instance().registerPair( std::string("UNIT"), CFG_TYPE_UNIT );
      success &= theCfgTypeMap::Instance().registerPair( std::string("CPPACS"), CFG_TYPE_CPPACS );
      success &= theCfgTypeMap::Instance().registerPair( std::string("CERN"), CFG_TYPE_CERN );
      success &= theCfgTypeMap::Instance().registerPair( std::string("WEAK_FIELD"), CFG_TYPE_WEAK_FIELD);
      success &= theCfgTypeMap::Instance().registerPair( std::string("CLASSICAL_SF"), CFG_TYPE_CLASSICAL_SF);
      
      return success;
    }

    bool registered = registerAll();
    const std::string typeIDString = "CfgType";
  }

  using namespace CfgTypeEnv;

  void read(XMLReader& xml_in, const std::string& path, CfgType& t) {
    theCfgTypeMap::Instance().read(typeIDString, xml_in,path,t);
  }
  
  void write(XMLWriter& xml_out, const std::string& path, const CfgType& t) {
    theCfgTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
}
