#ifndef enum_gaugeacttype_io_h
#define enum_gaugeacttype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"


using namespace std;
using namespace Chroma;

namespace Chroma {

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! GaugeAct type
  enum GaugeActType {
    GAUGE_ACT_TYPE_WILSON
  };


  namespace GaugeActTypeEnv { 
    extern const string typeIDString;
    extern const bool registered; 
    const bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<GaugeActType> > theGaugeActTypeMap;

  // Reader and writer

  //! Read a GaugeActType enum
  void read(XMLReader& r, const string& path, GaugeActType& t);

  //! Write an GaugeActType enum
  void write(XMLWriter& w, const string& path, const GaugeActType& t);

};
#endif
