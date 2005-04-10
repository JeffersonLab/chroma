// -*- C++ -*-
// $Id: wilson_gaugeact.h,v 1.5 2005-04-10 22:32:38 edwards Exp $
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

  /*! @ingroup gaugeacts */
  namespace WilsonGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct WilsonGaugeActParams 
  {
    // Base Constructor
    WilsonGaugeActParams();
    
    // Read params from some root path
    WilsonGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    AnisoParam_t aniso;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, WilsonGaugeActParams& param);
  

  //! Wilson gauge action
  /*! \ingroup gaugeacts
   *
   * The standard Wilson gauge action
   */

  class WilsonGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    WilsonGaugeAct(Handle< GaugeBC > gbc_, 
		   const Real& beta)
      {param.beta = beta; init(gbc_);}

    //! General GaugeBC
    WilsonGaugeAct(Handle< GaugeBC > gbc_, 
		   const Real& beta,
		   const AnisoParam_t& aniso)
      {param.beta = beta; param.aniso = aniso; init(gbc_);}

    //! Read beta from a param struct
    WilsonGaugeAct(Handle< GaugeBC > gbc_, 
		   const WilsonGaugeActParams& p) :
      param(p) {init(gbc_);}

    //! Copy constructor
    WilsonGaugeAct(const WilsonGaugeAct& a) : 
      param(a.param) {}

    //! Assignment
    WilsonGaugeAct& operator=(const WilsonGaugeAct& a)
    {param=a.param; return *this;}

    //! Is anisotropy used?
    bool anisoP() const {return param.aniso.anisoP;}

    //! Anisotropy factor
    const Real anisoFactor() const {return param.aniso.xi_0;}

    //! Anisotropic direction
    int tDir() const {return param.aniso.t_dir;}

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
    const Real getBeta(void) const {return param.beta;}

  protected:
    //! Private initializer
    void init(Handle< GaugeBC > gbc);

  private:
    Handle<PlaqGaugeAct> plaq;  // Hold a plaquette gaugeact
    WilsonGaugeActParams param; // parameters
  };

} // end namespace Chroma


#endif
