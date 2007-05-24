// -*- C++ -*-
// $Id: aniso_sym_spatial_gaugeact.h,v 3.2 2007-05-24 19:36:05 bjoo Exp $
/*! \file
 *  \brief Spatial part of Tree level LW action
 */

#ifndef __aniso_sym_spatial_gaugeact_h__
#define __aniso_sym_spatial_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_sym_gaugeact_params.h"


namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace AnisoSymSpatialGaugeActEnv 
  {
    extern const string name;
    bool registerAll();
  }

  //! Spatial anisotropic Symanzik improved gauge action
  /*! \ingroup gaugeacts
   *
   *   Contains space-space plaquette and space space rectangle terms
   *   only. Useful for when one wants to split the spatial and temporal
   *   parts of the general Symanzik gauge action onto different timescales
   *   in an (R)HMC simulation.
   */

  class AnisoSymSpatialGaugeAct : public LinearGaugeAction
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Read beta from a param struct
    AnisoSymSpatialGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
			  const AnisoSymSpatialGaugeActParams& p) :
      cgs(cgs_), param(p) {init();}

    //! Is anisotropy used?
    bool anisoP() const {return param.aniso.anisoP; }

    //! Anisotropy factor
    const Real anisoFactor() const {return param.aniso.xi_0;}

    //! Anisotropic direction
    int tDir() const {return param.aniso.t_dir;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const
    {
      QDPIO::cerr << "This function is not implemented" << endl;
      QDP_abort(1);
      
    }

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const;
      /*
	{
	plaq->derivSpatial(result,state,param.aniso.t_dir); 
	multi1d<LatticeColorMatrix> tmp;
	rect->derivSpatial(tmp,state);
	result += tmp;
	}
      */
    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const;
    /*
    {
      return plaq->spatialS(state, param.aniso.t_dir) + rect->spatialS(state);
    }
    */

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

    //! Destructor is automatic
    ~AnisoSymSpatialGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return param.beta; }
    
    const Real getUS(void) const { return param.u_s; }
    

  protected:
    //! Private initializer
    //    void init(Handle< CreateGaugeState<P,Q> > cgs_);
    void init(void);

    //! Hide assignment
    void operator=(const AnisoSymSpatialGaugeAct& a) {}

  private:
    const Handle< CreateGaugeState<P,Q> > cgs;
    AnisoSymSpatialGaugeActParams    param;    /*!< The couplings and anisotropy*/
    //   Handle<PlaqGaugeAct>           plaq;     /*!< Hold a plaquette gaugeact */
    //   Handle<RectGaugeAct>           rect;     /*!< Hold a rectangle gaugeact */

    Real plaq_c_s;
    Real rect_c_s;
    
  };

};


#endif
