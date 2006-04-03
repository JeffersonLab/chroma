// -*- C++ -*-
// $Id: schr_sf_gaugebc.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! @file
 * @brief Schroedinger gauge boundary conditions
 */

#ifndef __schr_sf_gaugebc_h__
#define __schr_sf_gaugebc_h__

#include "actions/gauge/gaugebcs/schroedinger_gaugebc.h"

namespace Chroma
{

  //! Abstract class for SOME Schroedinger gauge BC
  /*! @ingroup gaugebc
   *
   *  Schroedinger BC implies periodic in dirs orthog to decay dir, and some
   *  kind of fixed BC in the decay dir.
   */
  class SchrSFGaugeBC : public SchrGaugeBC
  {
  public:
    //! Virtual destructor
    virtual ~SchrSFGaugeBC() {}

    //! Decay direction
    virtual int getDir() const = 0;

  protected:
    //! Construct the mask and boundary fields
    virtual void initBnd(multi1d<LatticeColorMatrix>& SFBndFld,
			 multi1d<LatticeBoolean>& lSFmask) const;

    //! Maximum plaquette size. This is what knows about 1x1 plaq or 1x2 rect.
    /*! \return 1 for 1x1 plaq or 2 for 1x2 rect in decay_dir */
    virtual int getMaxExtent() const = 0;

    //! Multiplier on phases
    virtual const Real& SchrPhiMult() const = 0;

    //! Structure holding phases
    struct Phases_t
    {
      multi1d<Real> lower;
      multi1d<Real> upper;
    };

    //! Get the angles on the boundaries
    virtual const Phases_t& getPhases() const = 0;
  };

}



#endif
