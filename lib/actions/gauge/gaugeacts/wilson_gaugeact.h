// -*- C++ -*-
// $Id: wilson_gaugeact.h,v 1.1 2005-01-13 02:02:38 edwards Exp $
/*! \file
 *  \brief Wilson gauge action
 */

#ifndef __wilson_gaugeact_h__
#define __wilson_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

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
      gbc(gbc_), beta(beta_) {}

    //! Read beta from a param struct
    WilsonGaugeAct(Handle< GaugeBC > gbc_, 
		   const WilsonGaugeActParams& p) :
      gbc(gbc_), beta(p.beta) {}

    //! Copy constructor
    WilsonGaugeAct(const WilsonGaugeAct& a) : 
      gbc(a.gbc), beta(a.beta) {}

    //! Assignment
    WilsonGaugeAct& operator=(const WilsonGaugeAct& a)
    {gbc=a.gbc; beta=a.beta; return *this;}

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
    const GaugeBC& getGaugeBC() const {return *gbc;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		Handle<const ConnectState> state,
		int mu, int cb) const;

    //! Compute dS/dU
    void dsdu(multi1d<LatticeColorMatrix>& result,
	      const Handle<const ConnectState> state) const;

    //! Compute the actions
    Double S(const Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~WilsonGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return beta; }

  private:
    Handle< GaugeBC >  gbc;  // Gauge Boundary Condition
    Real beta;               // The coupling Beta

  };

};

using namespace Chroma;

#endif
