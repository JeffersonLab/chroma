// $Id: prec_constdet_two_flavor_monomial5d_w.cc,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/prec_constdet_two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"

#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      return new EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D(
	TwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }

    const std::string name("TWO_FLAVOR_EOPREC_CONSTDET_FERM_MONOMIAL5D");

    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;

      foo &= WilsonTypeFermActs5DEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(name, createMonomial);

      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  }; //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv



  // Constructor
  EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D::EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D(
    const TwoFlavorWilsonTypeFermMonomialParams& param_) 
  {
    inv_param = param_.inv_param;

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

    WilsonTypeFermAct5D<T,P,Q>* tmp_act = TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");

    EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>* downcast=dynamic_cast<EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecConstDetWilsonTypeFermAct5D in EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5D()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    // Get Chronological predictor
    AbsChronologicalPredictor5D<LatticeFermion>* tmp=0x0;
    if( param_.predictor_xml == "" ) {
      // No predictor specified use zero guess
      tmp = new ZeroGuess5DChronoPredictor(fermact->size());
    }
    else {
      try { 
	std::string chrono_name;
	std::istringstream chrono_is(param_.predictor_xml);
	XMLReader chrono_xml(chrono_is);
	read(chrono_xml, "/ChronologicalPredictor/Name", chrono_name);
	tmp = The5DChronologicalPredictorFactory::Instance().createObject(chrono_name, fermact->size(), chrono_xml, "/ChronologicalPredictor");
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
  }
  
} //end namespace Chroma
