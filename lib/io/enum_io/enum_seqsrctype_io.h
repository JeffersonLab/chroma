#ifndef enum_seqsrctype_io_h
#define enum_seqsrctype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"


using namespace std;
using namespace Chroma;

namespace Chroma {

  /*!
   * Types and structures
   *y
   * \ingroup io
   *
   * @{
   */
  //! Sequenti\al Source type
  enum SeqSourceType {
    SEQ_SOURCE_TYPE_NUCL_U_UNPOL  = 0,
    SEQ_SOURCE_TYPE_NUCL_D_UNPOL  = 1,
    SEQ_SOURCE_TYPE_NUCL_U_POL    = 2,   // same as polarized in 2pt baryon case
    SEQ_SOURCE_TYPE_NUCL_D_POL    = 3,   // same as polarized in 2pt baryon case
    SEQ_SOURCE_TYPE_DELTA_U_UNPOL = 4,
    SEQ_SOURCE_TYPE_DELTA_D_UNPOL = 5,
    SEQ_SOURCE_TYPE_NUCL_U_UNPOL_NONREL = 6,
    SEQ_SOURCE_TYPE_NUCL_D_UNPOL_NONREL = 7,
    SEQ_SOURCE_TYPE_NUCL_U_POL_NONREL   = 8,
    SEQ_SOURCE_TYPE_NUCL_D_POL_NONREL   = 9,
    SEQ_SOURCE_TYPE_NUCL_U_MIXED_NONREL = 21,
    SEQ_SOURCE_TYPE_NUCL_D_MIXED_NONREL = 22,
    SEQ_SOURCE_TYPE_PION = 10
  };


  namespace SeqSourceTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<SeqSourceType> > theSeqSourceTypeMap;

  // Reader and writer

  //! Read an SeqSourceType enum
  void read(XMLReader& r, const string& path, SeqSourceType& t);

  //! Write an SeqSourceType enum
  void write(XMLWriter& w, const string& path, const SeqSourceType& t);

};
#endif
