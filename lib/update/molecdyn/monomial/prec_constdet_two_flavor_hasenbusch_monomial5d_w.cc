// $Id: prec_constdet_two_flavor_hasenbusch_monomial5d_w.cc,v 3.2 2006-07-04 02:55:52 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_hasenbusch_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"

#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5DEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << name << endl;

      return new EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D(
	TwoFlavorHasenbuschWilsonTypeFermMonomialParams(xml, path));
    }
    
    const std::string name("TWO_FLAVOR_EOPREC_CONSTDET_HASENBUSCH_FERM_MONOMIAL5D");
 
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
  EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D::EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5D(
    const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param) 
  {
    inv_param = param.inv_param;

    std::istringstream is(param.fermact.xml);
    XMLReader fermact_reader(is);

    std::istringstream is_prec(param.fermact_prec.xml);
    XMLReader fermact_reader_prec(is_prec);
    
    if( param.fermact_prec.id != param.fermact.id ) {
      QDPIO::cerr << "For now both the numerator and the denominator fermacts mast be the same: You have asked for " 
		  << param.fermact.id 
		  << " in the denominator and " 
		  << param.fermact_prec.id << " in the numerator" << endl;
      QDP_abort(1);
      
    }

    QDPIO::cout << "EvenOddPrecTwoFlavorHasenbuschWilsonTypeFermMonomial5D: construct " 
		<< param.fermact.id << endl;
      
    WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
      TheWilsonTypeFermAct5DFactory::Instance().createObject(param.fermact.id, 
							     fermact_reader, 
							     param.fermact.path);

    EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>* downcast = 
      dynamic_cast<EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);
      
    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }
    
    fermact = downcast;    
      
    QDPIO::cout << "EvenOddPrecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " 
		<< param.fermact_prec.id << endl;
    
    WilsonTypeFermAct5D<T,P,Q>* tmp_act_prec = 
      TheWilsonTypeFermAct5DFactory::Instance().createObject(param.fermact_prec.id, 
							     fermact_reader_prec, 
							     param.fermact_prec.path);

    EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>* downcast_prec = 
      dynamic_cast<EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>*>(tmp_act_prec);
      
    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecConstDetWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }
    
    fermact_prec = downcast_prec;    
    
    if (fermact->size() != fermact_prec->size()) 
    {
      QDPIO::cerr << "Error: numerator action has to have the same length in the 5th dimension as the denominator action." << endl;
      QDPIO::cerr << "N5 in FermionAction " << fermact->size() << endl;
      QDPIO::cerr << "N5 in FermionActionPrec " << fermact_prec->size() << endl;
      QDP_abort(1);
    }
    
    // Get Chronological predictor
    AbsChronologicalPredictor5D<LatticeFermion>* tmp = 0x0;
    if ( param.predictor.xml == "" ) {
      // No predictor specified use zero guess
       tmp = new ZeroGuess5DChronoPredictor(fermact->size());
    }
    else 
    {
      try 
      { 
	std::istringstream chrono_is(param.predictor.xml);
	XMLReader chrono_xml(chrono_is);
	tmp = The5DChronologicalPredictorFactory::Instance().createObject(
	  param.predictor.id, fermact->size(),
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
  }

} //end namespace Chroma


