// -*- C++ -*-
// $Id: schr_nonpert_fermbc_w.h,v 3.2 2006-09-20 20:28:00 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - use for non-pertubative tuning of clover action
 */

#ifndef __schr_nonpert_fermbc_w_h__
#define __schr_nonpert_fermbc_w_h__

#include "actions/gauge/gaugebcs/schr_nonpert_gaugebc.h"
#include "actions/ferm/fermbcs/schr_sf_fermbc_w.h"
#include "actions/ferm/fermbcs/schr_fermbc_params_w.h"

namespace Chroma 
{ 

  /*! @ingroup fermbcs */
  namespace SchrNonPertFermBCEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Concrete class for Schroedinger BC - use for nonpertubative tuning
  /*! @ingroup fermbcs
   *
   *  Schroedinger BC for ferm actions
   */
  class SchrNonPertFermBC : public SchrSFFermBC
  {
  public:
    //! Only full constructor
    SchrNonPertFermBC(const SchrNonPertGaugeBC& gaugebc, 
		      const SchrFermBCParams& p) : param(p)
    {
      START_CODE();

      // Initialize boundary fields
      initBnd(fld, mask, maskF, gaugebc.SFBndFld(), gaugebc.lSFmask());

      END_CODE();
    }

    //! Destructor is automatic
    ~SchrNonPertFermBC() {}

    //! Decay direction
    int getDir() const {return param.decay_dir;}

  protected:
    //! Mask which lattice sites have fixed ferm sites
    const LatticeBoolean& lSFmaskF() const {return maskF;}

    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lSFmask() const {return mask;}

    //! Fixed gauge links on only the lSFmask() sites
    const multi1d<LatticeColorMatrix>& SFBndFld() const {return fld;}

    //! Maximum plaquette size. This is what knows about 1x1 plaq or 1x2 rect.
    /*! \return 1 for 1x1 plaq or 2 for 1x2 rect in decay_dir */
    int getMaxExtent() const {return param.loop_extent;}

    //! Get the angles on the boundaries
    const multi1d<Real>& getTheta() const {return param.theta;}

  private:
    // Hide default constuctor
    SchrNonPertFermBC() {}
    void operator=(const SchrNonPertFermBC&) {}

  private:
    SchrFermBCParams             param;
    LatticeBoolean               maskF;
    multi1d<LatticeBoolean>      mask;
    multi1d<LatticeColorMatrix>  fld;
  };
  
}

#endif
