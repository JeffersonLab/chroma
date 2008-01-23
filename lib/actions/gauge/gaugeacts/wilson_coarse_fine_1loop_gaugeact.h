// -*- C++ -*-
// $Id: wilson_coarse_fine_1loop_gaugeact.h,v 3.1 2008-01-23 15:36:11 edwards Exp $
/*! \file
 *  \brief Wilson gauge action supporting 2+2 anisotropy.
 *
 * Wilson gauge action on a 2+2 lattice.
 * Follows the conventions of  hep-lat/0303005 (TrinLat).
 * Here, coefficients are for 1-loop.
 */

#ifndef __wilson_coarse_fine_1loop_gaugeact_h__
#define __wilson_coarse_fine_1loop_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace WilsonCoarseFine1LoopGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct WilsonCoarseFine1LoopGaugeActParams 
  {
    //! Base Constructor
    WilsonCoarseFine1LoopGaugeActParams() {}
    
    //! Read params from some root path
    WilsonCoarseFine1LoopGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real  beta;      /*!< Yep, it's beta - the coefficient in front of the action */

    multi1d<int> coarse_dirs;  /*!< Directions for coarse lattice */
    Real  xi;                  /*!< Coefficient of fine-fine plaquette */
    Real  coeff_ff;            /*!< 1-loop coefficient of fine-fine plaquette = eta_ff-eta_cf */
    Real  coeff_cc;            /*!< 1-loop coefficient of coarse-coare plaquette = eta_cc-eta_cf */
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, WilsonCoarseFine1LoopGaugeActParams& param);
  

  //! WilsonCoarseFine gauge action
  /*! \ingroup gaugeacts
   *
   * Wilson gauge action on a 2+2 lattice.
   * Follows the conventions of  hep-lat/0303005 (TrinLat)
   */

  class WilsonCoarseFine1LoopGaugeAct : public LinearGaugeAction
  {
  public:
    //! Read beta from a param struct
    WilsonCoarseFine1LoopGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
				  const WilsonCoarseFine1LoopGaugeActParams& p);

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const
    {
      plaq->staple(result,state,mu,cb);
    }

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const
    {
      plaq->deriv(result,state);
    }

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const
    {
      return plaq->S(state);
    }

    //! Destructor is automatic
    ~WilsonCoarseFine1LoopGaugeAct() {}

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return plaq->getCreateState();}

  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! Partial constructor
    WilsonCoarseFine1LoopGaugeAct() {}
    //! Hide assignment
    void operator=(const WilsonCoarseFine1LoopGaugeAct& a) {}

  private:
    Handle<PlaqGaugeAct> plaq;             /*!< Hold a plaquette gaugeact */
    WilsonCoarseFine1LoopGaugeActParams param;  /*!< parameters */
  };

} // end namespace Chroma


#endif
