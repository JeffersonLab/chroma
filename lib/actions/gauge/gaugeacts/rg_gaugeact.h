// -*- C++ -*-
// $Id: rg_gaugeact.h,v 1.1 2005-01-13 04:30:51 edwards Exp $
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

  namespace RGGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  // Parameter structure
  struct RGGaugeActParams {
    // Base Constructor
    RGGaugeActParams();
    
    // Read params from some root path
    RGGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    Real c1;  
  };
  
  void read(XMLReader& xml, const string& path, RGGaugeActParams& param);
  

  //! RG gauge action
  /*! \ingroup gaugeact
   *
   * A RG (plaquette + rectangle) gauge action
   */

  class RGGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    RGGaugeAct(Handle< GaugeBC > gbc_, 
	       const Real& beta_, const Real& c1_) : 
      beta(beta_), c1(c1_) {init(gbc_);}

    //! Read beta from a param struct
    RGGaugeAct(Handle< GaugeBC > gbc_, 
	       const RGGaugeActParams& p) :
      beta(p.beta), c1(p.c1) {init(gbc_);}

    //! Copy constructor
    RGGaugeAct(const RGGaugeAct& a) : 
      beta(a.beta), c1(a.c1), plaq(a.plaq), rect(a.rect) {}

    //! Assignment
    RGGaugeAct& operator=(const RGGaugeAct& a)
    {beta=a.beta; c1=a.c1; plaq=a.plaq; rect=a.rect; return *this;}

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
    ~RGGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< GaugeBC > gbc);

  private:
    Real   beta;               // The coupling Beta
    Real   c1;                 // The relative rectangle factor
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact
    Handle<RectGaugeAct> rect; // Hold a rectangle gaugeact

  };

};

using namespace Chroma;

#endif
