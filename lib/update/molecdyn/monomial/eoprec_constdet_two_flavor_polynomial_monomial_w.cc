// $Id: eoprec_constdet_two_flavor_polynomial_monomial_w.cc,v 3.2 2008-05-23 21:31:33 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_polynomial_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomialEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
      {
	return new EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial(
	  TwoFlavorWilsonTypeFermMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    //! Identifier
    const std::string name = "TWO_FLAVOR_EOPREC_CONSTDET_POLYNOMIAL_FERM_MONOMIAL";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActs4DEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  } //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv


  // Constructor
  EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial::EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial(
    const TwoFlavorWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;

    {
      std::istringstream is(param.fermact.xml);
      XMLReader fermact_reader(is);
      
      QDPIO::cout << EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomialEnv::name 
		  << ": construct " << param.fermact.id << endl;
      
      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.fermact.id, fermact_reader, param.fermact.path);
      
      PolyWilsonTypeFermAct<T,P,Q>* downcast = 
	dynamic_cast<PolyWilsonTypeFermAct<T,P,Q>*>(tmp_act);
      
      // Check success of the downcast 
      if( downcast == 0x0 ) {
	QDPIO::cerr << "Unable to downcast FermAct to PolyWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial()" << endl;
	QDP_abort(1);
      }
      
      fermact = downcast;    
    }
    
    END_CODE();
  }

} //end namespace Chroma


