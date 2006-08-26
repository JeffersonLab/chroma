// $Id: prec_constdet_two_flavor_polyprec_monomial_w.cc,v 3.2 2006-08-26 02:08:42 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_polyprec_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"

#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomialEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << name << endl;

      return new EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial(
	TwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }
 
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;

      foo &= WilsonTypeFermActs4DEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(name, createMonomial);

      return foo;
    }

    //! Identifier
    const std::string name = "TWO_FLAVOR_EOPREC_CONSTDET_POLYPREC_FERM_MONOMIAL";

    //! Register the fermact
    const bool registered = registerAll();
  }; //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv


  // Constructor
  EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial::EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial(
    const TwoFlavorWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;

    {
      std::istringstream is(param.fermact.xml);
      XMLReader fermact_reader(is);
      
      QDPIO::cout << EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomialEnv::name 
		  << ": construct " << param.fermact.id << endl;
      
      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.fermact.id, fermact_reader, param.fermact.path);
      
      PolyWilsonTypeFermAct<T,P,Q>* downcast = 
	dynamic_cast<PolyWilsonTypeFermAct<T,P,Q>*>(tmp_act);
      
      // Check success of the downcast 
      if( downcast == 0x0 ) {
	QDPIO::cerr << "Unable to downcast FermAct to PolyWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomial()" << endl;
	QDP_abort(1);
      }
      
      fermact = downcast;    
    }

    // Get Chronological predictor
    AbsChronologicalPredictor4D<LatticeFermion>* tmp = 0x0;
    if( param.predictor.xml == "" ) {
      // No predictor specified use zero guess
       tmp = new ZeroGuess4DChronoPredictor();
    }
    else 
    {
      try 
      { 
	std::istringstream chrono_is(param.predictor.xml);
	XMLReader chrono_xml(chrono_is);
	tmp = The4DChronologicalPredictorFactory::Instance().createObject(param.predictor.id, 
									  chrono_xml, 
									  param.predictor.path);
      }
      catch(const std::string& e ) { 
	QDPIO::cerr << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }
    }
     
    if( tmp == 0x0 ) { 
      QDPIO::cerr << "Failed to create ZeroGuess4DChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;

    END_CODE();
  }
  
} //end namespace Chroma


