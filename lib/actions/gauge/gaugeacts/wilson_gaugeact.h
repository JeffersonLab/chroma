// -*- C++ -*-
// $Id: wilson_gaugeact.h,v 1.2 2005-01-13 04:29:45 edwards Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#ifndef __wilson_gaugeact_h__
#define __wilson_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"

namespace Chroma
{

  namespace WilsonGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  // Parameter structure
  struct WilsonGaugeActParams {
    // Base Constructor
    WilsonGaugeActParams();
    
    // Read params from some root path
    WilsonGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
  };
  
  void read(XMLReader& xml, const string& path, WilsonGaugeActParams& param);
  

  //! Wilson gauge action
  /*! \ingroup gaugeact
   *
   * The standard Wilson gauge action
   */

  class WilsonGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    WilsonGaugeAct(Handle< GaugeBC > gbc_, 
		   const Real& beta_) : 
      beta(beta_) {init(gbc_);}

    //! Read beta from a param struct
    WilsonGaugeAct(Handle< GaugeBC > gbc_, 
		   const WilsonGaugeActParams& p) :
      beta(p.beta) {init(gbc_);}

    //! Copy constructor
    WilsonGaugeAct(const WilsonGaugeAct& a) : 
      beta(a.beta), plaq(a.plaq) {}

    //! Assignment
    WilsonGaugeAct& operator=(const WilsonGaugeAct& a)
    {beta=a.beta; plaq=a.plaq; return *this;}

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
    }

    //! Compute dS/dU
    void dsdu(multi1d<LatticeColorMatrix>& result,
	      const Handle<const ConnectState> state) const
    {
      plaq->dsdu(result,state);
    }

    //! Compute the actions
    Double S(const Handle<const ConnectState> state) const
    {
      return plaq->S(state);
    }

    //! Destructor is automatic
    ~WilsonGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  protected:
    //! Private initializer
    void init(Handle< GaugeBC > gbc);

  private:
    Real   beta;               // The coupling Beta
    Handle<PlaqGaugeAct> plaq; // Hold a plaquette gaugeact

  };

};

using namespace Chroma;

#endif
