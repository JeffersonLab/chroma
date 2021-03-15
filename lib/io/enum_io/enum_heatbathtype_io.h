// -*- C++ -*-
/*! \file
 * \brief Heatbath enum
 */

#ifndef enum_heatbatbtype_io_h
#define enum_heatbatbtype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include "update/heatbath/su3hb.h"


namespace Chroma {

  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! Heatbath  type
  namespace HeatbathTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef Chroma::SingletonHolder<EnumTypeMap<HeatbathType> > theHeatbathTypeMap;

  // Reader and writer

  //! Read an HeatbathType enum
  void read(XMLReader& r, const std::string& path, HeatbathType& t);

  //! Write an HeatbathType enum
  void write(XMLWriter& w, const std::string& path, const HeatbathType& t);

  /*! @} */   // end of group io
}
#endif
