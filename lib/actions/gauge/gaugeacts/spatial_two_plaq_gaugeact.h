// -*- C++ -*-
// $Id: spatial_two_plaq_gaugeact.h,v 1.3 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#ifndef __spatial_two_plaq_gaugeact_h__
#define __spatial_two_plaq_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "io/aniso_io.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace SpatialTwoPlaqGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct SpatialTwoPlaqGaugeActParams 
  {
    // Base Constructor
    SpatialTwoPlaqGaugeActParams() {}
    
    // Read params from some root path
    SpatialTwoPlaqGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff;  
    AnisoParam_t aniso;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, SpatialTwoPlaqGaugeActParams& param);
  

  //! Plaquette gauge action
  /*! \ingroup gaugeacts
   *
   * The standard Plaquette gauge action
   */

  class SpatialTwoPlaqGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState<P,Q>
    SpatialTwoPlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const Real& coeff,
		 const AnisoParam_t& aniso) : 
      cgs(cgs_) {param.coeff = coeff; param.aniso = aniso; init();}

    //! Read coeff from a param struct
    SpatialTwoPlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const SpatialTwoPlaqGaugeActParams& p) :
      cgs(cgs_), param(p) {init();}

    //! Is anisotropy used?
    bool anisoP() const {return param.aniso.anisoP;}

    //! Anisotropy factor
    const Real anisoFactor() const {return param.aniso.xi_0;}

    //! Anisotropic direction
    int tDir() const {return param.aniso.t_dir;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
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
    ~SpatialTwoPlaqGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeff(void) const {return param.coeff;}

  protected:
    SpatialTwoPlaqGaugeAct() {}
    void operator=(const SpatialTwoPlaqGaugeAct& a) {}       //! Hide assignment

    //! Internal initializer
    void init();

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;  // Create Gauge State
    SpatialTwoPlaqGaugeActParams  param;   // The parameters

  };

};


#endif
