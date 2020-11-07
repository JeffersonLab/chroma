// -*- C++ -*-
/*! \file
 * \brief Wavetype enum
 */
#include "enum_wavetype_io.h"

#include <string>

namespace Chroma { 


  namespace WaveStateTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theWaveStateTypeMap::Instance().registerPair(std::string("S_WAVE"),
							   WAVE_TYPE_S_WAVE);

      success &=theWaveStateTypeMap::Instance().registerPair(std::string("P_WAVE" ), 
							  WAVE_TYPE_P_WAVE);
      
      success &=theWaveStateTypeMap::Instance().registerPair(std::string("D_WAVE" ), 
							  WAVE_TYPE_D_WAVE);

      return success;
    }

    const std::string typeIDString = "WaveStateType";
    bool registered = registerAll();
  }
  
  using namespace WaveStateTypeEnv;
  //! Read an WaveType enum
  void read(XMLReader& xml_in,  const std::string& path, WaveStateType& t) {
    theWaveStateTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an WaveType enum
  void write(XMLWriter& xml_out, const std::string& path, const WaveStateType& t) {
    theWaveStateTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
}
