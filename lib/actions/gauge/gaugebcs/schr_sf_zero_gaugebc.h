// -*- C++ -*-
// $Id: schr_sf_zero_gaugebc.h,v 3.1 2006-09-21 18:43:26 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - happens to zero out gauge fields in bc_dir
 */

#ifndef __schr_sf_zero_gaugebc_h__
#define __schr_sf_zero_gaugebc_h__

#include "actions/gauge/gaugebcs/schr_sf_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{ 

  /*! @ingroup gaugebcs */
  namespace SchrSFZeroGaugeBCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Concrete class for Schroedinger BC - zero out gauge boundaries
  /*! @ingroup gaugebcs
   *
   *  Schroedinger BC for gauge actions
   */
  class SchrSFZeroGaugeBC : public SchrSFGaugeBC
  {
  public:
    //! Only full constructor
    SchrSFZeroGaugeBC(const SchrGaugeBCParams& p);

    //! Destructor is automatic
    ~SchrSFZeroGaugeBC() {}

    //! Decay direction
    int getDir() const {return param.decay_dir;}

    //! Modify U fields in place
    /*! Default version provided */
    void modify(multi1d<LatticeColorMatrix>& u) const;

    //! Zero the some gauge-like field in place on the masked links
    /*! Default version provided */
    void zero(multi1d<LatticeColorMatrix>& ds_u) const;

    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lSFmask() const {return mask;}

    //! Fixed gauge links on only the lSFmask() sites
    const multi1d<LatticeColorMatrix>& SFBndFld() const {return fld;}

  protected:
    //! Maximum plaquette size. This is what knows about 1x1 plaq or 1x2 rect.
    /*! \return 1 for 1x1 plaq or 2 for 1x2 rect in decay_dir */
    int getMaxExtent() const {return param.loop_extent;}

    //! Multiplier on phases
    const Real& SchrPhiMult() const {return param.SchrPhiMult;}

    //! Get the angles on the boundaries
    const Phases_t& getPhases() const {return phases;}

    //! Initialize the phases
    void initPhases();

  private:
    // Hide default constuctor
    SchrSFZeroGaugeBC() {}
    void operator=(const SchrSFZeroGaugeBC&) {}

  private:
    SchrGaugeBCParams            param;
    Phases_t                     phases;
    multi1d<LatticeBoolean>      mask;
    multi1d<LatticeColorMatrix>  fld;
  };
  
}

#endif
