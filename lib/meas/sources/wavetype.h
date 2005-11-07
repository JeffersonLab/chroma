// -*- C++ -*-
// $Id: wavetype.h,v 2.1 2005-11-07 06:26:49 edwards Exp $
/*! \file
 *  \brief Particle wave type, like S-wave, P-Wave, etc.
 */

#ifndef __wavetype_h__
#define __wavetype_h__

namespace Chroma 
{ 

  /*! @ingroup sources */
  enum WaveStateType
  {
    WAVE_TYPE_S_WAVE, 
    WAVE_TYPE_P_WAVE, 
    WAVE_TYPE_D_WAVE
  };
}

#endif
