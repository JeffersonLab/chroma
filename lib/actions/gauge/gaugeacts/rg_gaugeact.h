// -*- C++ -*-
// $Id: rg_gaugeact.h,v 3.2 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief Generic RG style plaquette + rectangle gauge action
 */

#ifndef __rg_gaugeact_h__
#define __rg_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace RGGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct RGGaugeActParams {
    // Base Constructor
    RGGaugeActParams();
    
    // Read params from some root path
    RGGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    Real c1;  
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, RGGaugeActParams& param);
  

  //! RG gauge action
  /*! \ingroup gaugeacts
   *
   * A RG (plaquette + rectangle) gauge action
   */

  class RGGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState<P,Q>
    RGGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
	       const Real& beta_, const Real& c1_) : 
      beta(beta_), c1(c1_) {init(cgs_);}

    //! Read beta from a param struct
    RGGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
	       const RGGaugeActParams& p) :
      beta(p.beta), c1(p.c1) {init(cgs_);}

    //! Is anisotropy used?
    bool anisoP() const {return false;}

    //! Anisotropy factor
    const Real anisoFactor() const {return Real(1);}

    //! Anisotropic direction
    int tDir() const {return Nd-1;}

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

      LatticeColorMatrix tmp;
      rect->staple(tmp,state,mu,cb);
      result += tmp;
    }

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	      const Handle< GaugeState<P,Q> >& state) const
    {
      plaq->deriv(result,state);

      multi1d<LatticeColorMatrix> tmp;
      rect->deriv(tmp,state);
      result += tmp;
    }

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const
    {
      return plaq->S(state) + rect->S(state);
    }

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return plaq->getCreateState();}

    //! Destructor is automatic
    ~RGGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! Partial constructor
    RGGaugeAct() {}
    //! Hide assignment
    void operator=(const RGGaugeAct& a) {}

  private:
    Real   beta;               // The coupling Beta
    Real   c1;                 // The relative rectangle factor
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact
    Handle<RectGaugeAct> rect; // Hold a rectangle gaugeact

  };

};


#endif
