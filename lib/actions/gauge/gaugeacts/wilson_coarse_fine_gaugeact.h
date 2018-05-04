// -*- C++ -*-
/*! \file
 *  \brief Wilson gauge action supporting 2+2 anisotropy.
 *
 * Wilson gauge action on a 2+2 lattice.
 * Follows the conventions of  hep-lat/0303005 (TrinLat)
 */

#ifndef __wilson_coarse_fine_gaugeact_h__
#define __wilson_coarse_fine_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace WilsonCoarseFineGaugeActEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! @ingroup gaugeacts */
  struct WilsonCoarseFineGaugeActParams 
  {
    //! Base Constructor
    WilsonCoarseFineGaugeActParams() {}
    
    //! Read params from some root path
    WilsonCoarseFineGaugeActParams(XMLReader& xml_in, const std::string& path);

    multi1d<int> coarse_dirs;  /*!< Directions for coarse lattice */
    Real  coeff_ff;            /*!< General coefficient of fine-fine plaquette */
    Real  coeff_cf;            /*!< General coefficient of coarse-fine plaquette */
    Real  coeff_cc;            /*!< General coefficient of coarse-coare plaquette */
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const std::string& path, WilsonCoarseFineGaugeActParams& param);
  

  //! WilsonCoarseFine gauge action
  /*! \ingroup gaugeacts
   *
   * Wilson gauge action on a 2+2 lattice.
   * Follows the conventions of  hep-lat/0303005 (TrinLat)
   */

  class WilsonCoarseFineGaugeAct : public LinearGaugeAction
  {
  public:
    //! Read beta from a param struct
    WilsonCoarseFineGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
			     const WilsonCoarseFineGaugeActParams& p);

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
    ~WilsonCoarseFineGaugeAct() {}

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return plaq->getCreateState();}

  protected:
    //! Private initializer
    void init(Handle< CreateGaugeState<P,Q> > cgs);

    //! Partial constructor
    WilsonCoarseFineGaugeAct() {}
    //! Hide assignment
    void operator=(const WilsonCoarseFineGaugeAct& a) {}

  private:
    Handle<PlaqGaugeAct> plaq;             /*!< Hold a plaquette gaugeact */
    WilsonCoarseFineGaugeActParams param;  /*!< parameters */
  };

} // end namespace Chroma


#endif
