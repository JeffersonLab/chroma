// $Id: unprec_two_flavor_hasenbusch_monomial5d_w.cc,v 3.1 2006-06-21 20:42:09 bjoo Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5DEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      return new UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D(
	TwoFlavorHasenbuschWilsonTypeFermMonomialParams(xml, path));
    }
 
    const std::string name("TWO_FLAVOR_UNPREC_HASENBUSCH_FERM_MONOMIAL5D");
 
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
  } //end namespace Unprec TwoFlavorHasenbuschWilsonFermMonomialEnv


  // Constructor
  UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D::UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5D(
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
    XMLReader fermact_prec_reader(is_prec);

    // Get the name of the ferm act
    std::string fermact_prec_string;
    
    try { 
      read(fermact_prec_reader, "/FermionActionPrec/FermAct", fermact_prec_string);
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the fermact name: " << e<<  endl;
      QDP_abort(1);
    }

    // Check that the two 
    if( fermact_prec_string != fermact_string ) { 
      QDPIO::cerr << "For now both the numerator and the denominator fermacts mast be the same: You have asked for " 
		  << fermact_string 
		  << " in the denominator and " 
		  << fermact_prec_string << " in the numerator" << endl;
      QDP_abort(1);
    }
    

    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << fermact_string << endl;
    WilsonTypeFermAct5D<T,P,Q>* tmp_act = TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");
   

    UnprecWilsonTypeFermAct5D<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    


    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << fermact_prec_string << endl;

    WilsonTypeFermAct5D<T,P,Q>* tmp_act_prec = 
      TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact_prec_string, 
							     fermact_prec_reader, 
							     "/FermionActionPrec");

    UnprecWilsonTypeFermAct5D<T,P,Q>* downcast_prec = 
      dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act_prec);


    // Check success of the downcast 
    if( downcast_prec == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact_prec = downcast_prec;    

    if (fermact->size() != fermact_prec->size()) { 
      QDPIO::cerr << "Error: numerator action has to have the same length in the 5th dimension as the denominator action." << endl;
      QDPIO::cerr << "N5 in FermionAction " << fermact->size() << endl;
      QDPIO::cerr << "N5 in FermionActionPrec " << fermact_prec->size() << endl;
      QDP_abort(1);
    }
      

    //------------------------------------

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
	tmp = The5DChronologicalPredictorFactory::Instance().createObject(chrono_name, 
									  fermact->size(),
									  chrono_xml, 
									  std::string("/ChronologicalPredictor"));
      }
      catch(const std::string& e ) { 
	QDPIO::cerr << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }


    }
     
    if( tmp == 0x0 ) { 
      QDPIO::cerr << "Failed to create the 4D ChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;

    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial: finished " << fermact_string << endl;
  }
}; //end namespace Chroma


