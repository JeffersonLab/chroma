// $Id: unprec_two_flavor_monomial5d_w.cc,v 3.3 2006-09-20 20:28:05 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 5D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"

namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
      {
	return new UnprecTwoFlavorWilsonTypeFermMonomial5D(
	  TwoFlavorWilsonTypeFermMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name("TWO_FLAVOR_UNPREC_FERM_MONOMIAL5D");

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActs5DEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  } //end namespace Unprec TwoFlavorWilsonFermMonomialEnv



  // Constructor
  UnprecTwoFlavorWilsonTypeFermMonomial5D::UnprecTwoFlavorWilsonTypeFermMonomial5D(
    const TwoFlavorWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;

    std::istringstream is(param.fermact.xml);
    XMLReader fermact_reader(is);
    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial5D: construct " << param.fermact.id << endl;

    WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
      TheWilsonTypeFermAct5DFactory::Instance().createObject(param.fermact.id, fermact_reader, param.fermact.path);

    UnprecWilsonTypeFermAct5D<T,P,Q>* downcast = 
      dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct5D in UnprecTwoFlavorWilsonTypeFermMonomial5D()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    // Get Chronological predictor
    AbsChronologicalPredictor5D<LatticeFermion>* tmp = 0x0;
    if( param.predictor.xml == "" ) {
      // No predictor specified use zero guess
       tmp = new ZeroGuess5DChronoPredictor(fermact->size());
    }
    else 
    {
      try 
      { 
	std::string chrono_name;
	std::istringstream chrono_is(param.predictor.xml);
	XMLReader chrono_xml(chrono_is);
	tmp = The5DChronologicalPredictorFactory::Instance().createObject(
	  param.predictor.id, fermact->size(), chrono_xml, param.predictor.path);
      }
      catch(const std::string& e ) { 
	QDPIO::cerr << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }
    }
    
    if( tmp == 0x0 ) { 
      QDPIO::cerr << "Failed to create the 5D ChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;

    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial5D: finished " << param.fermact.id << endl;
    
    END_CODE();
  }

} //end namespace Chroma


