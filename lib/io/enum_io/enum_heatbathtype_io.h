#ifndef enum_heatbatbtype_io_h
#define enum_heatbatbtype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include "update/heatbath/su3hb.h"


namespace Chroma {

  namespace HeatbathTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<HeatbathType> > theHeatbathTypeMap;

  // Reader and writer

  //! Read an HeatbathType enum
  void read(XMLReader& r, const string& path, HeatbathType& t);

  //! Write an HeatbathType enum
  void write(XMLWriter& w, const string& path, const HeatbathType& t);

};
#endif
