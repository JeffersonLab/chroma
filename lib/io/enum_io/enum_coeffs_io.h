#ifndef enum_coeffs_io_h
#define enum_coeffs_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"


using namespace std;
using namespace Chroma;

namespace Chroma {
  // CfgType --------------------------------------
 
  enum CoeffType {
      COEFF_TYPE_ZOLOTAREV = 0,
      COEFF_TYPE_TANH
  };


  namespace CoeffTypeEnv { 
    extern const bool registered; 
    const bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<CoeffType> > theCoeffTypeMap;

  // Reader and writer
  void read(XMLReader& r, const string& path, CoeffType& t);
  void write(XMLWriter& w, const string& path, const CoeffType& t);

};
#endif
