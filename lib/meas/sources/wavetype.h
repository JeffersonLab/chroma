// -*- C++ -*-
// $Id: wavetype.h,v 3.0 2006-04-03 04:59:06 edwards Exp $
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
