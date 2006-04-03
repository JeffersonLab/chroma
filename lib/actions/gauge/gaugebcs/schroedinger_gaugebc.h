// -*- C++ -*-
// $Id: schroedinger_gaugebc.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
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
  class SchrGaugeBC : public GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! Virtual destructor
    virtual ~SchrGaugeBC() {}

    //! Modify U fields in place
    /*! Default version provided */
    void modify(multi1d<LatticeColorMatrix>& u) const;

    //! Zero the some gauge-like field in place on the masked links
    /*! Default version provided */
    void zero(multi1d<LatticeColorMatrix>& ds_u) const;

    //! Says if there are fixed links within the lattice
    bool nontrivialP() const {return true;}

    //! Decay direction
    virtual int getDir() const = 0;

    //! Mask which lattice sites have fixed gauge links
    virtual const multi1d<LatticeBoolean>& lSFmask() const = 0;

    //! Fixed gauge links on only the lSFmask() sites
    virtual const multi1d<LatticeColorMatrix>& SFBndFld() const = 0;
  };

}



#endif
