// -*- C++ -*-
// $Id: wavetype.h,v 1.2 2004-09-22 17:25:01 bjoo Exp $
/*! \file
 *  \brief Particle wave type, like S-wave, P-Wave, etc.
 */

#ifndef __wavetype_h__
#define __wavetype_h__

namespace Chroma { 

  enum WaveStateType
    {
      WAVE_TYPE_S_WAVE, 
      WAVE_TYPE_P_WAVE, 
      WAVE_TYPE_D_WAVE
    };
};

using namespace Chroma;
#endif
