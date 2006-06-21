// $Id: prec_constdet_two_flavor_hasenbusch_monomial5d_w.cc,v 3.1 2006-06-21 20:42:09 bjoo Exp $
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
    const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param_) 
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


    std::istringstream is_prec(param_.ferm_act_prec);
    XMLReader fermact_reader_prec(is_prec);
    
    // Get the name of the ferm act
    std::string fermact_string_prec;
    try { 
      read(fermact_reader_prec, "/FermionActionPrec/FermAct", fermact_string_prec);
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the fermact name: " << e<<  endl;
      QDP_abort(1);
    }
    
    if( fermact_string_prec != fermact_string ) {
      QDPIO::cerr << "For now both the numerator and the denominator fermacts mast be the same: You have asked for " 
		  << fermact_string 
		  << " in the denominator and " 
		  << fermact_string_prec << " in the numerator" << endl;
      QDP_abort(1);
      
    }

    QDPIO::cout << "EvanOddPrecTwoFlavorHasenbuschWilsonTypeFermMonomial5D: construct " << fermact_string << endl;
    
      
    WilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* tmp_act = TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");

    EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>* downcast=dynamic_cast<EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);
      
    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }
    
    fermact = downcast;    
    
      
    QDPIO::cout << "EvanOddPrecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << fermact_string_prec << endl;
    
    WilsonTypeFermAct5D<T,P,Q>* tmp_act_prec = TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionActionPrec");

    EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>* downcast_prec=dynamic_cast<EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>*>(tmp_act_prec);
      
    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecConstDetWilsonTypeFermAct in EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }
    
    fermact_prec = downcast_prec;    
    
    if (fermact->size() != fermact_prec->size()) { 
      QDPIO::cerr << "Error: numerator action has to have the same length in the 5th dimension as the denominator action." << endl;
      QDPIO::cerr << "N5 in FermionAction " << fermact->size() << endl;
      QDPIO::cerr << "N5 in FermionActionPrec " << fermact_prec->size() << endl;
      QDP_abort(1);
    }
    
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
	tmp = The5DChronologicalPredictorFactory::Instance().createObject(chrono_name, fermact->size(),
								 chrono_xml, 
								 "/ChronologicalPredictor");
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


  
}; //end namespace Chroma


