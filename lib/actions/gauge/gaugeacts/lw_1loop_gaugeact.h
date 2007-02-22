// -*- C++ -*-
// $Id: lw_1loop_gaugeact.h,v 3.2 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief 1-loop tadpole-improved Luscher-Weisz gauge action
 */

#ifndef __lw_1loop_gaugeact_h__
#define __lw_1loop_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/pg_gaugeact.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace LW1LoopGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct LW1LoopGaugeActParams {
    // Base Constructor
    LW1LoopGaugeActParams();
    
    // Read params from some root path
    LW1LoopGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    Real u0;  
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, LW1LoopGaugeActParams& param);
  

  //! LW1Loop gauge action
  /*! \ingroup gaugeacts
   *
   * The standard LW1Loop gauge action
   */

  class LW1LoopGaugeAct : public LinearGaugeAction
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General GaugeState
    LW1LoopGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		    const Real& beta_, const Real& u0_) : 
      beta(beta_), u0(u0_) {init(cgs_);}

    //! Read beta from a param struct
    LW1LoopGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		    const LW1LoopGaugeActParams& p) :
      beta(p.beta), u0(p.u0) {init(cgs_);}

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

      pg->staple(tmp,state,mu,cb);
      result += tmp;
    }

    //! Compute dS/dU
    void deriv(P& result,
	       const Handle< GaugeState<P,Q> >& state) const
    {
      plaq->deriv(result,state);

      multi1d<LatticeColorMatrix> tmp;
      rect->deriv(tmp,state);
      result += tmp;

      pg->deriv(tmp,state);
      result += tmp;
    }

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const
    {
      return plaq->S(state) + rect->S(state) + pg->S(state);
    }

    //! Produce a gauge boundary condition object
    const CreateGaugeState<P,Q>& getCreateState() const {return plaq->getCreateState();}

    //! Destructor is automatic
    ~LW1LoopGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! General GaugeState
    LW1LoopGaugeAct() {}
    //! Hide assignment
    void operator=(const LW1LoopGaugeAct& a) {}

  private:
    Real   beta;               // The coupling Beta
    Real   u0;                 // Tadpole factor
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact
    Handle<RectGaugeAct> rect; // Hold a rectangle gaugeact
    Handle<PgGaugeAct>   pg;   // Hold a parallelogram gaugeact

  };

};


#endif
