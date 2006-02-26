// -*- C++ -*-
// $Id: schr2link_gaugebc.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Schroedinger functional 2-link gauge boundary conditions
 */

#ifndef __gaugebc_schr2link_h__
#define __gaugebc_schr2link_h__

#include "gaugebc.h"
#include "actions/gauge/gaugebcs/schroedinger_gaugebc.h"
#include "actions/gauge/gaugebcs/schr_gaugebc_params.h"

namespace Chroma 
{ 

  /*! @ingroup gaugebcs */
  namespace Schr2LinkGaugeBCEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Concrete class for 2-link gauge action boundary conditions with Schroedinger BC
  /*! @ingroup gaugebcs
   *
   *  Schroedinger BC for gauge actions that have only 2 links for padding
   *  in the decay direction - e.g., Symanzik-like gauge action
   */
  class Schr2LinkGaugeBC : public SchrGaugeBC
  {
  public:
    //! Schroedinger BC
    /*! 
     * \param SchrFun_      type of Schroedinger BC
     * \param SchrPhiMult_  factor to rescale fixed field
     */
    Schr2LinkGaugeBC(const SchrFunType SchrFun_, const Real& SchrPhiMult_)
      {
	initFunc(SchrFun_, SchrPhiMult_);
      }

    //! Initialise from params struct
    Schr2LinkGaugeBC(const SchrGaugeBCParams& params)
      {
	initFunc(params.SchrFun, params.SchrPhiMult);
      }

    //! Copy constructor
    Schr2LinkGaugeBC(const Schr2LinkGaugeBC& a) : SchrFun(a.SchrFun), 
      decay_dir(a.decay_dir), mask(a.mask), fld(a.fld) {}

    //! Destructor is automatic
    ~Schr2LinkGaugeBC() {}

    //! Assignment
    Schr2LinkGaugeBC& operator=(const Schr2LinkGaugeBC& a)
      {
	SchrFun = a.SchrFun; decay_dir = a.decay_dir; 
	mask = a.mask; fld = a.fld;
	return *this;
      }

#if defined(EXPOSE_THIS_STUFF)
    //! Type of Schroedinger BC
    SchrFunType getSFBC() const {return SchrFun;}
#endif

    //! Modify U fields in place
    void modify(multi1d<LatticeColorMatrix>& u) const
      {QDP_error_exit("modify not implemented");}

    //! Zero the U fields in place on the masked links
    void zero(multi1d<LatticeColorMatrix>& u) const
      {QDP_error_exit("zero not implemented");}

#if defined(EXPOSE_THIS_STUFF)
    // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lbmaskU() const {return mask;}

    //! Fixed gauge links on only the lbmaskU() sites
    const multi1d<LatticeColorMatrix>& lFldU() const {return fld;}
#endif

    //! Says if there are fixed links within the lattice
    bool nontrivialP() const {return true;}

    //! Decay direction
    int getDir() const {return decay_dir;}

  private:
    // Hide default constuctor
    Schr2LinkGaugeBC() {}

    void initFunc(const SchrFunType SchrFun_, const Real& SchrPhiMult_) {
      QDPIO::cerr << "Schr2LinkGaugeBC() not yet implemented" << endl;
      QDP_abort(1);
    }

  private:
    SchrFunType SchrFun;
    int decay_dir;
    multi1d<LatticeBoolean> mask;
    multi1d<LatticeColorMatrix> fld;
  };

  
}

#endif
