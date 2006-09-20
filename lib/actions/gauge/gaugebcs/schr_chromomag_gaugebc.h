// -*- C++ -*-
// $Id: schr_chromomag_gaugebc.h,v 3.1 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - chromo-magnetic gauge BC
 */

#ifndef __schr_chromomag_gaugebc_h__
#define __schr_chromomag_gaugebc_h__

#include "actions/gauge/gaugebcs/schroedinger_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{ 

  /*! @ingroup gaugebcs */
  namespace SchrChromoMagGaugeBCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Concrete class for Schroedinger BC - use for nonpertubative tuning
  /*! @ingroup gaugebcs
   *
   *  Schroedinger BC for gauge actions
   */
  class SchrChromoMagGaugeBC : public SchrGaugeBC
  {
  public:
    //! Only full constructor
    SchrChromoMagGaugeBC(const SchrGaugeBCParams& p);

    //! Destructor is automatic
    ~SchrChromoMagGaugeBC() {}

    //! Decay direction
    int getDir() const {return param.decay_dir;}

    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lSFmask() const {return mask;}

    //! Fixed gauge links on only the lSFmask() sites
    const multi1d<LatticeColorMatrix>& SFBndFld() const {return fld;}

  private:
    // Hide default constuctor
    SchrChromoMagGaugeBC() {}
    void operator=(const SchrChromoMagGaugeBC&) {}

  private:
    SchrGaugeBCParams            param;
    multi1d<LatticeBoolean>      mask;
    multi1d<LatticeColorMatrix>  fld;
  };
  
}

#endif
