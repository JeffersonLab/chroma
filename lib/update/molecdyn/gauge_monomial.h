// -*- C++ -*-
// $Id: gauge_monomial.h,v 1.1 2005-01-12 04:43:33 edwards Exp $
/*! \file
 *  \brief Generic gauge action monomial wrapper
 */

#ifndef __gaugeact_monomial_h__
#define __gaugeact_monomial_h__

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
#include "actions/gauge/gaugeact_factory.h"

namespace Chroma 
{
  namespace GaugeMonomialEnv {
    extern const string name;
    extern const bool registered;
  };

  // Parameter structure
  struct GaugeMonomialParams 
  {
    // Base Constructor
    GaugeMonomialParams();

    // Read monomial from some root path
    GaugeMonomialParams(XMLReader& in, const std::string& path);
    string gauge_act;
  };


  //! Wrapper class for  gauge monomials
  /*!
   * Monomial is expected to be the same for these fermacts
   */
  class GaugeMonomial :
    public  ExactMonomial< multi1d<LatticeColorMatrix>,
                           multi1d<LatticeColorMatrix> >    
    {
    public: 
      // Construct out of a parameter struct. Check against the desired GaugeAct name
      GaugeMonomial(const string& gaugeact_name, 
		    const GaugeMonomialParams& param_);

      // Copy Constructor
      GaugeMonomial(const GaugeMonomial& m) : gaugeact((m.gaugeact)) {}

      // Create a suitable state and compute F
      void dsdq(multi1d<LatticeColorMatrix>& F, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) {
	// Make a gauge connect state
	Handle< const ConnectState> g_state(getGaugeAct().createState(s.getQ()));

	getGaugeAct().dsdu(F, g_state);
      }


      Double S(const AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) const {
	Handle< const ConnectState> g_state(getGaugeAct().createState(s.getQ()));
	return getGaugeAct().S(g_state);
      }
	
	
      void refreshInternalFields(const AbsFieldState<multi1d<LatticeColorMatrix>,
				 multi1d<LatticeColorMatrix> >& s) {
	//No internal fields to refresh => Nop
      }

      void setInternalFields(const Monomial<multi1d<LatticeColorMatrix>, 
			     multi1d<LatticeColorMatrix> >& m) {
	// No internal fields to refresh => Nop
      }

    protected:
      const GaugeAction& getGaugeAct(void) const { 
	return *gaugeact;
      }

    private:
      // Hide empty constructor and =
      GaugeMonomial();
      void operator=(const GaugeMonomial&);

    private:
      // A handle for the UnprecFermAct
      Handle<const GaugeAction> gaugeact;
    };


}; //end namespace chroma













#endif
