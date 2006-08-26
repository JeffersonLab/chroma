// $Id: unprec_two_flavor_hasenbusch_monomial_w.cc,v 3.2 2006-08-26 02:08:42 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorHasenbuschWilsonTypeFermMonomialEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      return new UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(
	TwoFlavorHasenbuschWilsonTypeFermMonomialParams(xml, path));
    }
 
    const std::string name("TWO_FLAVOR_UNPREC_HASENBUSCH_FERM_MONOMIAL");
 
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;

      foo &= WilsonTypeFermActs4DEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(name, createMonomial);

      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  } //end namespace Unprec TwoFlavorHasenbuschWilsonFermMonomialEnv


  // Constructor
  UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial::UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(
    const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;

    std::istringstream is(param.fermact.xml);
    XMLReader fermact_reader(is);
    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << param.fermact.id << endl;

   
    WilsonTypeFermAct<T,P,Q>* tmp_act = 
      TheWilsonTypeFermActFactory::Instance().createObject(param.fermact.id, fermact_reader, param.fermact.path);

    UnprecWilsonTypeFermAct<T,P,Q>* downcast = 
      dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);


    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    // -------------------
    std::istringstream is_prec(param.fermact_prec.xml);
    XMLReader fermact_prec_reader(is_prec);
    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << param.fermact_prec.id << endl;

    WilsonTypeFermAct<T,P,Q>* tmp_act_prec = 
      TheWilsonTypeFermActFactory::Instance().createObject(param.fermact_prec.id, 
							   fermact_prec_reader, 
							   param.fermact_prec.path);

    UnprecWilsonTypeFermAct<T,P,Q>* downcast_prec = 
      dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act_prec);


    // Check success of the downcast 
    if( downcast_prec == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact_prec = downcast_prec;    

    //------------------------------------

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
      QDPIO::cerr << "Failed to create the 4D ChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;

    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial: finished " << param.fermact.id << endl;
    
    END_CODE();
  }

} //end namespace Chroma


