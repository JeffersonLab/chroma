#include "chromabase.h"

#include "update/molecdyn/wilson_gauge_monomial.h"
#include "update/molecdyn/monomial_factory.h"
#include "actions/gauge/wilson_gaugeact.h"

#include <string>

using namespace QDP;
using namespace std;

namespace Chroma { 
 
  namespace WilsonGaugeMonomialEnv {
    //! Callback function for the factory
    ExactMonomial< multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> >*
		      createMonomial(XMLReader& xml, const string& path) {

      XMLReader gaugebc_xml(xml, path);
      QDPIO::cout << "Create Monomial: " << WilsonGaugeMonomialEnv::name << endl;

      GaugeAction *gaugeact =  TheGaugeActFactory::Instance().createObject(WilsonGaugeActEnv::name, gaugebc_xml, "./GaugeAction");

      // Downcast
      const WilsonGaugeAct* w_gaugeact=dynamic_cast<const WilsonGaugeAct*>(gaugeact);
      if( w_gaugeact == 0 ) { 
	QDPIO::cerr << "Unable to create WilsonGaugeAct " << endl;
	QDP_abort(1);
      }
      Handle< const WilsonGaugeAct > w_handle(w_gaugeact);

      return new WilsonGaugeMonomial(w_handle);
    }
    
    const std::string name="WILSON_GAUGE_MONOMIAL";
    const bool registered=TheExactMonomialFactory::Instance().registerObject(name, createMonomial);
    
  }; //end namespace Unprec TwoFlavorWilsonFermMonomialEnv




}; //end namespace Chroma


