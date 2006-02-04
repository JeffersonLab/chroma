// -*- C++ -*-
// $Id: enum_wavetype_io.h,v 2.1 2006-02-04 16:42:25 edwards Exp $
/*! \file
 * \brief Wavetype enum
 */
#ifndef enum_wavetype_io_h
#define enum_wavetype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"
#include "meas/sources/wavetype.h"


namespace Chroma {

  /*************** SOURCES ********************************/
  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */
  //! Wave state type

  namespace WaveStateTypeEnv { 
    extern const string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<WaveStateType> > theWaveStateTypeMap;

  //! Read an WaveStateType enum
  void read(XMLReader& r, const string& path, WaveStateType& t);

  //! Write an WaveStateType enum
  void write(XMLWriter& w, const string& path, const WaveStateType& t);

  /*! @} */   // end of group io

};
#endif
