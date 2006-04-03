// $Id: prec_constdet_two_flavor_polynomial_monomial_w.cc,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_polynomial_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomialEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << name << endl;

      return new EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial(
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
    const std::string name = "TWO_FLAVOR_EOPREC_CONSTDET_POLYNOMIAL_FERM_MONOMIAL";

    //! Register the fermact
    const bool registered = registerAll();
  }; //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv


  // Constructor
  EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial::EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial(
    const TwoFlavorWilsonTypeFermMonomialParams& param_) 
  {
    inv_param = param_.inv_param;

    {
      std::istringstream is(param_.ferm_act);
      XMLReader fermact_reader(is);
      
      // Get the name of the ferm act
      std::string fermact_string;
      try { 
	read(fermact_reader, "/FermionAction/FermAct", fermact_string);
      }
      catch( const std::string& e) { 
	QDPIO::cerr << "Error grepping the fermact name: " << e<<  endl;
	QDP_abort(1);
      }
      
      QDPIO::cout << EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomialEnv::name << ": construct " << fermact_string << endl;
      
      WilsonTypeFermAct<T,P,Q>* tmp_act = TheWilsonTypeFermActFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");
      
      PolyWilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<PolyWilsonTypeFermAct<T,P,Q>*>(tmp_act);
      
      // Check success of the downcast 
      if( downcast == 0x0 ) {
	QDPIO::cerr << "Unable to downcast FermAct to PolyWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomial()" << endl;
	QDP_abort(1);
      }
      
      fermact = downcast;    
    }
    
  }


  
}; //end namespace Chroma


