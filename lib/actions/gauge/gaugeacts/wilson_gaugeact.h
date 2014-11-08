// -*- C++ -*-
// $Id: wilson_gaugeact.h,v 3.4 2007-03-22 19:06:26 bjoo Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#ifndef __wilson_gaugeact_h__
#define __wilson_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/wilson_gaugeact_params.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace WilsonGaugeActEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

  //! Wilson gauge action
  /*! \ingroup gaugeacts
   *
   * The standard Wilson gauge action
   */

  class WilsonGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState<P,Q>
    WilsonGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		   const Real& beta)
      {param.beta = beta; init(cgs_);}

    //! General CreateGaugeState<P,Q>
    WilsonGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		   const Real& beta,
		   const AnisoParam_t& aniso)
      {param.beta = beta; param.aniso = aniso; init(cgs_);}

    //! Read beta from a param struct
    WilsonGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		   const WilsonGaugeActParams& p) :
      param(p) {init(cgs_);}

    //! Is anisotropy used?
    bool anisoP() const {return param.aniso.anisoP;}

    //! Anisotropy factor
    const Real anisoFactor() const {return param.aniso.xi_0;}

    //! Anisotropic direction
    int tDir() const {return param.aniso.t_dir;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const
    {
      plaq->staple(result,state,mu,cb);
    }

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const
    {
      plaq->deriv(result,state);
    }

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const
    {
      return plaq->S(state);
    }

    //! Destructor is automatic
    ~WilsonGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const {return param.beta;}

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return plaq->getCreateState();}

  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! Partial constructor
    WilsonGaugeAct() {}
    //! Hide assignment
    void operator=(const WilsonGaugeAct& a) {}

  private:
    Handle<PlaqGaugeAct> plaq;  // Hold a plaquette gaugeact
    WilsonGaugeActParams param; // parameters
  };

} // end namespace Chroma


#endif
