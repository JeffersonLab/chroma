#include "enum_wavetype_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 


  namespace WaveStateTypeEnv { 

    const bool registerAll(void) 
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

    const bool registered = registerAll();
  };
  
  //! Read an WaveType enum
  void read(XMLReader& xml_in,  const string& path, WaveStateType& t) {
    theWaveStateTypeMap::Instance().read(xml_in, path,t);
  }
  
  //! Write an WaveType enum
  void write(XMLWriter& xml_out, const string& path, const WaveStateType& t) {
    theWaveStateTypeMap::Instance().write(xml_out, path, t);
  }
};
