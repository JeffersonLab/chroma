// -*- C++ -*-
// $Id: wavetype.h,v 2.0 2005-09-25 21:04:40 edwards Exp $
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

#endif
