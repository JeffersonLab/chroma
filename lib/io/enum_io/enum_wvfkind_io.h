#ifndef enum_wvfkind_io_h
#define enum_wvfkinf_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include "meas/smear/wvfkind.h"  // This is where the enum is defined

using namespace std;
using namespace Chroma;

namespace Chroma {

  namespace WvfKindEnv { 
    extern const string typeIDString;
    extern const bool registered; 
    const bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<WvfKind> > theWvfKindMap;

  // Reader and writer

  //! Read an WvfKind enum
  void read(XMLReader& r, const string& path, WvfKind& t);

  //! Write an WvfKind Type enum
  void write(XMLWriter& w, const string& path, const WvfKind& t);

};
#endif
