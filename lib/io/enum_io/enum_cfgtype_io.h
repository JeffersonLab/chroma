#ifndef enum_cfgtype_io_h
#define enum_cfgtype_io_h


#include "io/enum_io/enum_type_map.h"
#include "singleton.h"
#include <string>

using namespace std;
using namespace Chroma;

namespace Chroma {
 
  //! Configuration type
  enum CfgType {
      CFG_TYPE_MILC = 0,
      CFG_TYPE_NERSC,
      CFG_TYPE_SCIDAC,
      CFG_TYPE_SZIN,
      CFG_TYPE_SZINQIO,
      CFG_TYPE_DISORDERED,
      CFG_TYPE_UNIT
  };


  namespace CfgTypeEnv { 
    extern const bool registered; 
    const bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<CfgType> > theCfgTypeMap;

  //! read a configuration type enum
  void read(XMLReader& xml_in, const string& path, CfgType& t);

  //! write a configuration type enum
  void write(XMLWriter& xml_out, const string& path, const CfgType& t); 

};

#endif
