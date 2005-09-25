// -*- C++ -*-
// $Id: lw_tree_gaugeact.h,v 2.0 2005-09-25 21:04:31 edwards Exp $
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

  class LWTreeGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    LWTreeGaugeAct(Handle< GaugeBC > gbc_, 
		   const Real& beta_, const Real& u0_) : 
      beta(beta_), u0(u0_) {init(gbc_);}

    //! Read beta from a param struct
    LWTreeGaugeAct(Handle< GaugeBC > gbc_, 
		   const LWTreeGaugeActParams& p) :
      beta(p.beta), u0(p.u0) {init(gbc_);}

    //! Copy constructor
    LWTreeGaugeAct(const LWTreeGaugeAct& a) : 
      beta(a.beta), u0(a.u0), plaq(a.plaq), rect(a.rect) {}

    //! Assignment
    LWTreeGaugeAct& operator=(const LWTreeGaugeAct& a)
    {beta=a.beta; u0=a.u0; plaq=a.plaq; rect=a.rect; return *this;}

    //! Is anisotropy used?
    bool anisoP() const {return false;}

    //! Anisotropy factor
    const Real anisoFactor() const {return Real(1);}

    //! Anisotropic direction
    int tDir() const {return Nd-1;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const OrderedSet& getSet() const {return rb;}

    //! Produce a gauge boundary condition object
    const GaugeBC& getGaugeBC() const {return plaq->getGaugeBC();}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		Handle<const ConnectState> state,
		int mu, int cb) const
    {
      plaq->staple(result,state,mu,cb);

      LatticeColorMatrix tmp;
      rect->staple(tmp,state,mu,cb);
      result += tmp;
    }

    //! Compute dS/dU
    void dsdu(multi1d<LatticeColorMatrix>& result,
	      const Handle<const ConnectState> state) const
    {
      plaq->dsdu(result,state);

      multi1d<LatticeColorMatrix> tmp;
      rect->dsdu(tmp,state);
      result += tmp;
    }

    //! Compute the actions
    Double S(const Handle<const ConnectState> state) const
    {
      return plaq->S(state) + rect->S(state);
    }

    //! Destructor is automatic
    ~LWTreeGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< GaugeBC > gbc);

  private:
    Real   beta;               // The coupling Beta
    Real   u0;                 // Tadpole factor
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact
    Handle<RectGaugeAct> rect; // Hold a rectangle gaugeact

  };

};


#endif
