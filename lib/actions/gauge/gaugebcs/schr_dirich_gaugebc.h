// -*- C++ -*-
// $Id: schr_dirich_gaugebc.h,v 3.1 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - dirichlet gauge BC
 */

#ifndef __schr_dirich_gaugebc_h__
#define __schr_dirich_gaugebc_h__

#include "actions/gauge/gaugebcs/schroedinger_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{ 

  /*! @ingroup gaugebcs */
  namespace SchrDirichletGaugeBCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Concrete class for Schroedinger BC - use for nonpertubative tuning
  /*! @ingroup gaugebcs
   *
   *  Schroedinger BC for gauge actions
   */
  class SchrDirichletGaugeBC : public SchrGaugeBC
  {
  public:
    //! Only full constructor
    SchrDirichletGaugeBC(const SchrGaugeBCParams& p);

    //! Destructor is automatic
    ~SchrDirichletGaugeBC() {}

    //! Decay direction
    int getDir() const {return param.decay_dir;}

    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lSFmask() const {return mask;}

    //! Fixed gauge links on only the lSFmask() sites
    const multi1d<LatticeColorMatrix>& SFBndFld() const {return fld;}

  private:
    // Hide default constuctor
    SchrDirichletGaugeBC() {}
    void operator=(const SchrDirichletGaugeBC&) {}

  private:
    SchrGaugeBCParams            param;
    multi1d<LatticeBoolean>      mask;
    multi1d<LatticeColorMatrix>  fld;
  };
  
}

#endif
