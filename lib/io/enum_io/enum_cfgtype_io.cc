#include "enum_cfgtype_io.h"
#include <string>

using namespace std;
using namespace Chroma;

namespace Chroma { 

  namespace CfgTypeEnv { 

    const bool registerAll(void) 
    {
      bool success; 
      success = theCfgTypeMap::Instance().registerPair(string("MILC"), CFG_TYPE_MILC );
      success &= theCfgTypeMap::Instance().registerPair(string("NERSC"), CFG_TYPE_NERSC);
      success &= theCfgTypeMap::Instance().registerPair( string("MILC"), CFG_TYPE_MILC );
      success &= theCfgTypeMap::Instance().registerPair( string("NERSC"), CFG_TYPE_NERSC );
      success &= theCfgTypeMap::Instance().registerPair( string("SCIDAC"), CFG_TYPE_SCIDAC );
      success &= theCfgTypeMap::Instance().registerPair( string("SZIN" ), CFG_TYPE_SZIN );
      success &= theCfgTypeMap::Instance().registerPair( string("SZINQIO"), CFG_TYPE_SZINQIO );
      success &= theCfgTypeMap::Instance().registerPair( string("DISORDERED"), CFG_TYPE_DISORDERED );
      success &= theCfgTypeMap::Instance().registerPair( string("UNIT"), CFG_TYPE_UNIT );
      
      return success;
    }

    const bool registered = registerAll();
  };

  void read(XMLReader& xml_in, const string& path, CfgType& t) {
    theCfgTypeMap::Instance().read(xml_in,path,t);
  }
  
  void write(XMLWriter& xml_out, const string& path, const CfgType& t) {
    theCfgTypeMap::Instance().write(xml_out, path, t);
  }
};
