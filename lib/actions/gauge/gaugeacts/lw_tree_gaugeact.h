// -*- C++ -*-
// $Id: lw_tree_gaugeact.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Tree-level tadpole-improved Luscher-Weisz gauge action
 */

#ifndef __lw_tree_gaugeact_h__
#define __lw_tree_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace LWTreeGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct LWTreeGaugeActParams {
    // Base Constructor
    LWTreeGaugeActParams();
    
    // Read params from some root path
    LWTreeGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    Real u0;  
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, LWTreeGaugeActParams& param);
  

  //! LWTree gauge action
  /*! \ingroup gaugeacts
   *
   * The standard LWTree gauge action
   */

  class LWTreeGaugeAct : public LinearGaugeAction
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General CreateGaugeState
    LWTreeGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		   const Real& beta_, const Real& u0_) : 
      beta(beta_), u0(u0_) {init(cgs_);}

    //! Read beta from a param struct
    LWTreeGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		   const LWTreeGaugeActParams& p) :
      beta(p.beta), u0(p.u0) {init(cgs_);}

    //! Is anisotropy used?
    bool anisoP() const {return false;}

    //! Anisotropy factor
    const Real anisoFactor() const {return Real(1);}

    //! Anisotropic direction
    int tDir() const {return Nd-1;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const OrderedSet& getSet() const {return rb;}

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
    ~LWTreeGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! Hide assignment
    void operator=(const LWTreeGaugeAct& a) {}

  private:
    Real   beta;               // The coupling Beta
    Real   u0;                 // Tadpole factor
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact
    Handle<RectGaugeAct> rect; // Hold a rectangle gaugeact

  };

};


#endif
