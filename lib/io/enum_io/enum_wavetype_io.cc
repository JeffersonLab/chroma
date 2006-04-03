// -*- C++ -*-
// $Id: enum_wavetype_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
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
      success = theWaveStateTypeMap::Instance().registerPair(string("S_WAVE"),
							   WAVE_TYPE_S_WAVE);

      success &=theWaveStateTypeMap::Instance().registerPair(string("P_WAVE" ), 
							  WAVE_TYPE_P_WAVE);
      
      success &=theWaveStateTypeMap::Instance().registerPair(string("D_WAVE" ), 
							  WAVE_TYPE_D_WAVE);

      return success;
    }

    const string typeIDString = "WaveStateType";
    bool registered = registerAll();
  };
  
  using namespace WaveStateTypeEnv;
  //! Read an WaveType enum
  void read(XMLReader& xml_in,  const string& path, WaveStateType& t) {
    theWaveStateTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an WaveType enum
  void write(XMLWriter& xml_out, const string& path, const WaveStateType& t) {
    theWaveStateTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
