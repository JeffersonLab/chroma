// -*- C++ -*-
// $Id: gaugebc.h,v 1.9 2004-12-30 10:29:36 bjoo Exp $
/*! @file
 * @brief Gauge boundary conditions
 */

#ifndef __gaugebc_h__
#define __gaugebc_h__

#include "chromabase.h"
#include "io/enum_io/enum_gaugebc_io.h"

using namespace QDP;

namespace Chroma
{
  //! Base class for all gauge action boundary conditions
  /*! @ingroup actions
   *
   */
  class GaugeBC
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~GaugeBC() {}

    //! Apply the BC onto the U fields in place
    virtual void modify(multi1d<LatticeColorMatrix>& u) const = 0;

    //! Zero the U fields in place on the masked links
    virtual void zero(multi1d<LatticeColorMatrix>& u) const = 0;

#if defined(EXPOSE_THIS_STUFF)
    // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
    //! Mask which lattice sites have fixed gauge links
    virtual const multi1d<LatticeBoolean>& lbmaskU() const = 0;

    //! Fixed gauge links on only the lbmaskU() sites
    virtual const multi1d<LatticeColorMatrix>& lFldU() const = 0;
#endif

    //! Says if there are fixed links within the lattice
    virtual bool nontrivialP() const = 0;
  };

  //! Abstract class for all gauge action boundary conditions with Schroedinger BC
  /*! @ingroup actions
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




 
};


using namespace Chroma;

#endif
