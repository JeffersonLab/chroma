// -*- C++ -*-
// $Id: schr_coupling_gaugebc.h,v 3.1 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - use for coupling determinations
 */

#ifndef __schr_coupling_gaugebc_h__
#define __schr_coupling_gaugebc_h__

#include "actions/gauge/gaugebcs/schr_sf_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{ 

  /*! @ingroup gaugebcs */
  namespace SchrCouplingGaugeBCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Concrete class for Schroedinger BC - use for coupling determination
  /*! @ingroup gaugebcs
   *
   *  Schroedinger BC for gauge actions 
   */
  class SchrCouplingGaugeBC : public SchrSFGaugeBC
  {
  public:
    //! Only full constructor
    SchrCouplingGaugeBC(const SchrGaugeBCParams& p);

    //! Destructor is automatic
    ~SchrCouplingGaugeBC() {}

    //! Decay direction
    int getDir() const {return param.decay_dir;}

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
    SchrCouplingGaugeBC() {}
    void operator=(const SchrCouplingGaugeBC&) {}

  private:
    SchrGaugeBCParams            param;
    Phases_t                     phases;
    multi1d<LatticeBoolean>      mask;
    multi1d<LatticeColorMatrix>  fld;
  };
  
}

#endif
