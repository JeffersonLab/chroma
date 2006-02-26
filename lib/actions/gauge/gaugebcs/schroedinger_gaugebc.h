// -*- C++ -*-
// $Id: schroedinger_gaugebc.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! @file
 * @brief Schroedinger Gauge boundary conditions
 */

#ifndef __schroedinger_gaugebc_h__
#define __schroedinger_gaugebc_h__

#include "gaugebc.h"

namespace Chroma
{

  //! Abstract class for all gauge action boundary conditions with Schroedinger BC
  /*! @ingroup gaugebc
   *
   *  Schroedinger BC implies periodic in dirs orthog to decay dir, and some
   *  kind of fixed BC in the decay dir.
   */
  class SchrGaugeBC : public GaugeBC
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SchrGaugeBC() {}

#if defined(EXPOSE_THIS_STUFF)
    //! Type of Schroedinger BC
    virtual SchrFunType getSFBC() const = 0;
#endif

    //! Decay direction
    virtual int getDir() const = 0;
  };

}



#endif
