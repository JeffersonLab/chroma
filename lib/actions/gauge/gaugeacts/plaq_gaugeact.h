// -*- C++ -*-
// $Id: plaq_gaugeact.h,v 1.3 2005-02-23 19:26:12 edwards Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#ifndef __plaq_gaugeact_h__
#define __plaq_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "io/aniso_io.h"

namespace Chroma
{

  namespace PlaqGaugeActEnv 
  { 
    extern const string name;
    extern const bool registered;
  }

  // Parameter structure
  struct PlaqGaugeActParams 
  {
    // Base Constructor
    PlaqGaugeActParams() {}
    
    // Read params from some root path
    PlaqGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff;  
    AnisoParam_t aniso;
  };
  
  void read(XMLReader& xml, const string& path, PlaqGaugeActParams& param);
  

  //! Plaquette gauge action
  /*! \ingroup gaugeact
   *
   * The standard Plaquette gauge action
   */

  class PlaqGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    PlaqGaugeAct(Handle< GaugeBC > gbc_, 
		 const Real& coeff,
		 const AnisoParam_t& aniso) : 
      gbc(gbc_) {param.coeff = coeff; param.aniso = aniso;}

    //! Read coeff from a param struct
    PlaqGaugeAct(Handle< GaugeBC > gbc_, 
		 const PlaqGaugeActParams& p) :
      gbc(gbc_), param(p) {}

    //! Copy constructor
    PlaqGaugeAct(const PlaqGaugeAct& a) : 
      gbc(a.gbc), param(a.param) {}

    //! Assignment
    PlaqGaugeAct& operator=(const PlaqGaugeAct& a)
    {gbc=a.gbc; param=a.param; return *this;}

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
    ~PlaqGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeff(void) const {return param.coeff;}

  private:
    Handle< GaugeBC >  gbc;  // Gauge Boundary Condition
    PlaqGaugeActParams  param; // The parameters

  };

};


#endif
