// -*- C++ -*-
// $Id: plaq_plus_spatial_two_plaq_gaugeact.h,v 3.5 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief Plaquette gauge action
 */

#ifndef __plaq_plus_spatial_two_plaq_gaugeact_h__
#define __plaq_plus_spatial_two_plaq_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "io/aniso_io.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace PlaqPlusSpatialTwoPlaqGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
    extern  double getTime();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct PlaqPlusSpatialTwoPlaqGaugeActParams 
  {
    // Base Constructor
    PlaqPlusSpatialTwoPlaqGaugeActParams() {}
    
    // Read params from some root path
    PlaqPlusSpatialTwoPlaqGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff_plaq_s;
    Real coeff_plaq_t;

    Real coeff_two_plaq;  
    AnisoParam_t aniso;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, PlaqPlusSpatialTwoPlaqGaugeActParams& param);
  

  //! Plaquette gauge action
  /*! \ingroup gaugeacts
   *
   * The standard Plaquette gauge action
   */

  class PlaqPlusSpatialTwoPlaqGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState<P,Q>
    PlaqPlusSpatialTwoPlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
				   const Real& coeff_plaq_s_,
				   const Real& coeff_plaq_t_,
				   const Real& coeff_two_plaq_,
		 const AnisoParam_t& aniso) : 
      cgs(cgs_) {
      param.coeff_plaq_s = coeff_plaq_s_; 
      param.coeff_plaq_t = coeff_plaq_t_;
      param.coeff_two_plaq = coeff_two_plaq_;
      param.aniso = aniso; 
      init();
    }

    //! Read coeff from a param struct
    PlaqPlusSpatialTwoPlaqGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		 const PlaqPlusSpatialTwoPlaqGaugeActParams& p) :
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
    ~PlaqPlusSpatialTwoPlaqGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeffPlaqS(void) const {return param.coeff_plaq_s;}
    const Real getCoeffPlaqT(void) const {return param.coeff_plaq_t;}
    const Real getCoeffTwoPlaq(void) const {return param.coeff_two_plaq;}

  protected:
    PlaqPlusSpatialTwoPlaqGaugeAct() {}
    void operator=(const PlaqPlusSpatialTwoPlaqGaugeAct& a) {}       //! Hide assignment

    //! Internal initializer
    void init();

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;  // Create Gauge State
    PlaqPlusSpatialTwoPlaqGaugeActParams  param;   // The parameters

  };

};


#endif
