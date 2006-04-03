// -*- C++ -*-
// $Id: enum_srcsnktype_io.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief SrcSink enum
 */
#ifndef enum_srcsnktype_io_h
#define enum_srcsnktype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include "meas/sources/srcsnktype.h"


namespace Chroma {

  /*************** SOURCES ********************************/
  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! Source  type
  namespace SourceTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<SourceType> > theSourceTypeMap;

  //! Read an SourceType enum
  void read(XMLReader& r, const string& path, SourceType& t);

  //! Write an SourceType enum
  void write(XMLWriter& w, const string& path, const SourceType& t);

  /*************  SINKS **********************************/

  namespace SinkTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<SinkType> > theSinkTypeMap;

  // Reader and writer

  //! Read an SinkType enum
  void read(XMLReader& r, const string& path, SinkType& t);

  //! Write an SinkType enum
  void write(XMLWriter& w, const string& path, const SinkType& t);

  /*! @} */   // end of group io

};
#endif
