// -*- C++ -*-
// $Id: mp_gaugeact.h,v 1.1 2006-07-20 18:03:17 bjoo Exp $
/*! \file
 *  \brief Morningstar Peardont type gaugeact
 */

#ifndef __mp_tree_gaugeact_h__
#define __mp_tree_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/two_plaq_spatial_gaugeact.h"


namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace MPTGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct MPGaugeActParams {
    // Base Constructor
    MPGaugeActParams();
    
    // Read params from some root path
    MPGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real beta;  
    Real u0s;
    Real u0t;
    Real omega;
    AnisoParam_t aniso;
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, MPGaugeActParams& param);
  

  //! MP gauge action
  /*! \ingroup gaugeacts
   *
   * The standard MP gauge action
   */

  class MPGaugeAct : public LinearGaugeAction
  {
  public:
    // Typedefs to save typing
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Read beta from a param struct
    MPGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
		   const MPGaugeActParams& p) :
      param(p) {init(cgs_);}

    //! Is anisotropy used?
    bool anisoP() const {return param.aniso.anisoP; }

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
      rect->staple(tmp,state,mu,cb);     // This may fail I think in aniso mode
      result += tmp;
      two_plaq->staple(tmp,state,mu,cb); // This will fail I think
      result += tmp;
    }

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	      const Handle< GaugeState<P,Q> >& state) const
    {
      plaq->deriv(result,state);

      multi1d<LatticeColorMatrix> tmp;
      rect->deriv(tmp,state);
      result += tmp;

      two_plaq->deriv(tmp,state);
      result += tmp;

    }

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const
    {
      return plaq->S(state) + rect->S(state) + two_plaq->S(state);
    }

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return plaq->getCreateState();}

    //! Destructor is automatic
    ~MPGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getBeta(void) const { return param.beta; }
    
    const Real getU0S(void) const { return param.u0s; }
    
    const Real getU0T(void) const { return param.u0t; }

    const Real getOmega(void) const { return param.omega; }


  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! Hide assignment
    void operator=(const MPGaugeAct& a) {}

  private:
    Handle<PlaqGaugeAct> plaq;    /*!< Hold a plaquette gaugeact */
    Handle<RectGaugeAct> rect;    /*!< Hold a rectangle gaugeact */
    Handle<TwoPlaqSpatialGaugeAct> two_plaq /*!< Hold spatial plaquettes separated in time type gaugeact */
    MPGaugeActParams param;   /*!< The couplings and anisotropy*/
  };

};


#endif
