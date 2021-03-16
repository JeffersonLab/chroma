// -*- C++ -*-
/*! \file
 * \brief Wavetype enum
 */
#ifndef enum_wavetype_io_h
#define enum_wavetype_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"


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
  enum WaveStateType
  {
    WAVE_TYPE_S_WAVE, 
    WAVE_TYPE_P_WAVE, 
    WAVE_TYPE_D_WAVE
  };

  namespace WaveStateTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered; 
    bool registerAll(void);   // Forward declaration
  }

  // A singleton to hold the typemap
  typedef Chroma::SingletonHolder<EnumTypeMap<WaveStateType> > theWaveStateTypeMap;

  //! Read an WaveStateType enum
  void read(XMLReader& r, const std::string& path, WaveStateType& t);

  //! Write an WaveStateType enum
  void write(XMLWriter& w, const std::string& path, const WaveStateType& t);

  /*! @} */   // end of group io

}
#endif
