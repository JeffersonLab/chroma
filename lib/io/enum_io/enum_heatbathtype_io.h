// -*- C++ -*-
// $Id: enum_heatbathtype_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
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

  /*! @} */   // end of group io
};
#endif
