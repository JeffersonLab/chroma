// -*- C++ -*-
// $Id: plaq_gaugeact.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
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

  /*! @ingroup gaugeacts */
  namespace PlaqGaugeActEnv 
  { 
    extern const string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct PlaqGaugeActParams 
  {
    // Base Constructor
    PlaqGaugeActParams() {}
    
    // Read params from some root path
    PlaqGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff;  
    AnisoParam_t aniso;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, PlaqGaugeActParams& param);
  

  //! Plaquette gauge action
  /*! \ingroup gaugeacts
   *
   * The standard Plaquette gauge action
   */

  class PlaqGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState<P,Q>
    PlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const Real& coeff,
		 const AnisoParam_t& aniso) : 
      cgs(cgs_) {param.coeff = coeff; param.aniso = aniso;}

    //! Read coeff from a param struct
    PlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const PlaqGaugeActParams& p) :
      cgs(cgs_), param(p) {}

    //! Is anisotropy used?
    bool anisoP() const {return param.aniso.anisoP;}

    //! Anisotropy factor
    const Real anisoFactor() const {return param.aniso.xi_0;}

    //! Anisotropic direction
    int tDir() const {return param.aniso.t_dir;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const OrderedSet& getSet() const {return rb;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const;

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const;

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const;

    //! Destructor is automatic
    ~PlaqGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeff(void) const {return param.coeff;}

  protected:
    PlaqGaugeAct() {}
    void operator=(const PlaqGaugeAct& a) {}       //! Hide assignment

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;  // Create Gauge State
    PlaqGaugeActParams  param; // The parameters

  };

};


#endif
