#ifndef enum_qdpvolfmt_io_h
#define enum_qdpvolfmt_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"


using namespace std;
using namespace Chroma;

namespace Chroma {


  namespace QDPVolfmtEnv { 
    extern const bool registered; 
    const bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<QDP_volfmt_t> > theQDPVolfmtMap;

  // Reader and writer

  //! Read a QDP volume format type
  void read(XMLReader& r, const string& path, QDP_volfmt_t& t);
  
  //! Write a QDP volume format type
  void write(XMLWriter& w, const string& path, const QDP_volfmt_t& t);

};
#endif
