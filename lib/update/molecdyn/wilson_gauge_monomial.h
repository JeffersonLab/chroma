#ifndef WILSON_GAUGE_MONOMIAL_H
#define WILSON_GAUGE_MONOMIAL_H

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
#include "actions/gauge/wilson_gaugeact.h"
#include "actions/gauge/gaugeact_factory.h"

using namespace std;
using namespace QDP;
using namespace std;

namespace Chroma {

  namespace WilsonGaugeMonomialEnv {
    extern const string name;
    extern const bool registered;
  };


  class WilsonGaugeMonomial :
    public  ExactMonomial< multi1d<LatticeColorMatrix>,
                           multi1d<LatticeColorMatrix> >    
    {
    public: 

      // Construct from a gaugeact handle
      WilsonGaugeMonomial(Handle< const WilsonGaugeAct > gaugeact_) : gaugeact(gaugeact_) {}

      // Copy Constructor
      WilsonGaugeMonomial(const WilsonGaugeMonomial& m) : gaugeact((m.gaugeact)) {}


      // Create a suitable state and compute F
      void dsdq(multi1d<LatticeColorMatrix>& F, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const {

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
      const WilsonGaugeAct& getGaugeAct(void) const { 
	return *gaugeact;
      }

    private:
 
      // A handle for the UnprecWilsonFermAct
      Handle<const WilsonGaugeAct> gaugeact;
    };


}; //end namespace chroma













#endif
