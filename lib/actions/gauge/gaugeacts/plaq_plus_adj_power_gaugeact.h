// -*- C++ -*-
/*! \file
 *  \brief Plaquette plus a power of an adjoint gauge action
 *
 *  Implements:
 *    S =  beta*Re(1 - (1/Nc)*Tr(U_p)) + gamma*Re((1 - (1/Nc^2)*Tr(U_p^dag)*Tr(U_p))^(q/2))
 */

#ifndef __plaq_plus_adj_power_gaugeact_h__
#define __plaq_plus_adj_power_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace PlaqPlusAdjPowerGaugeActEnv 
  { 
    bool registerAll();

    //! Parameter structure
    /*! @ingroup gaugeacts */
    struct Params 
    {
      //! Base Constructor
      Params() {}
    
      //! Read params from some root path
      Params(XMLReader& xml_in, const std::string& path);

      Real  beta;   // Coupling for fundamental plaquette
      Real  gamma;  // Coupling for adjoint power
      int   q;      // Power of plaquette
    };
  

    //! PlaqPlusAdjPower  gauge action
    /*! \ingroup gaugeacts
     *
     * The standard  gauge action as sum of PlaqPlusAdjPowers
     */
    class GaugeAct : public LinearGaugeAction
    {
    public:
      // Typedefs to save typing
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //! General CreateGaugeState<P,Q>
      //! Read coeff from a param struct
      GaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, const Params& p) :
	cgs(cgs_), param(p) {}

      //! Return the set on which the gauge action is defined
      /*! Defined on the even-off (red/black) set */
      const Set& getSet() const {return rb;}

      //! Compute staple
      /*! Default version. Derived class should override this if needed. */
      void staple(LatticeColorMatrix& result,
		  const Handle< GaugeState<P,Q> >& state,
		  int mu, int cb) const;

      //! Compute dS/dU
      void deriv(multi1d<LatticeColorMatrix>& result,
		 const Handle< GaugeState<P,Q> >& state) const;

      //! Compute the actions
      Double S(const Handle< GaugeState<P,Q> >& state) const;

      //! Produce a gauge create state object
      const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

      //! Destructor is automatic
      ~GaugeAct() {}

    protected:
      //! Hide assignment
      void operator=(const GaugeAct& a) {}
      
      //! Compute the site-level action
      void siteAction(multi2d<LatticeComplex>& site_act, 
		      const Handle< GaugeState<P,Q> >& state) const;

    private:
      Handle< CreateGaugeState<P,Q> >  cgs;  /*!< Create Gauge State */
      Params               param;            /*!< The parameters */
    };

  }

}


#endif
