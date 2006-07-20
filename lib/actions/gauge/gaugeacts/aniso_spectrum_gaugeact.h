// -*- C++ -*-
// $Id: aniso_spectrum_gaugeact.h,v 1.1 2006-07-20 18:40:35 edwards Exp $
/*! \file
 *  \brief Anisotropic gaugeact useful for spectrum from hep-lat/9911003
 *
 *  Tree-level LW with tapole improvement, missing 1x2 in time, also including
 *  2-plaq term. Taken from Morningstar-Peardon, hep-lat/9911003
 */

#ifndef __aniso_spectrum_gaugeact_h__
#define __aniso_spectrum_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/spatial_two_plaq_gaugeact.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace AnisoSpectrumGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct AnisoSpectrumGaugeActParams {
    // Base Constructor
    AnisoSpectrumGaugeActParams();
    
    // Read params from some root path
    AnisoSpectrumGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;             /*!< Coupling */
    Real u_s;              /*!< Spatial tadpole */
    Real u_t;              /*!< Temporal tadpole */
    Real omega;            /*!< Weighting of adjoint term */
    AnisoParam_t aniso;    /*!< Anisotropy struct */
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, AnisoSpectrumGaugeActParams& param);
  

  //! AnisoSpectrum gauge action
  /*! \ingroup gaugeacts
   *
   * The standard AnisoSpectrum gauge action
   */

  class AnisoSpectrumGaugeAct : public LinearGaugeAction
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Read beta from a param struct
    AnisoSpectrumGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		    const AnisoSpectrumGaugeActParams& p) :
      param(p) {init(cgs_);}

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
		int mu, int cb) const
    {
      plaq->staple(result,state,mu,cb);

      LatticeColorMatrix tmp;
      rect->staple(tmp,state,mu,cb);
      result += tmp;

      plaq_sq->staple(tmp,state,mu,cb);
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

      plaq_sq->deriv(tmp,state);
      result += tmp;
    }

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const
    {
      return plaq->S(state) + rect->S(state) + plaq_sq->S(state);
    }

    //! Produce a gauge boundary condition object
    const CreateGaugeState<P,Q>& getCreateState() const {return plaq->getCreateState();}

    //! Destructor is automatic
    ~AnisoSpectrumGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return param.beta; }

  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! General GaugeState
    AnisoSpectrumGaugeAct() {}

    //! Hide assignment
    void operator=(const AnisoSpectrumGaugeAct& a) {}

  private:
    AnisoSpectrumGaugeActParams    param;    /*!< Parameters */
    Handle<PlaqGaugeAct>           plaq;     /*!< Hold a plaquette gaugeact */
    Handle<RectGaugeAct>           rect;     /*!< Hold a rectangle gaugeact */
    Handle<SpatialTwoPlaqGaugeAct> plaq_sq;  /*!< Hold a spatial 2-plaq gaugeact */

  };

};


#endif
