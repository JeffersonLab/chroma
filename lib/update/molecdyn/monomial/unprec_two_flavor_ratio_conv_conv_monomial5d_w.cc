// $Id: unprec_two_flavor_ratio_conv_conv_monomial5d_w.cc,v 3.1 2008-05-23 21:31:36 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_ratio_conv_conv_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5DEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
      {
	return new UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5D(
	  TwoFlavorRatioConvConvWilsonTypeFermMonomialParams(xml, path));
      }
 
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("TWO_FLAVOR_UNPREC_RATIO_CONV_CONV_FERM_MONOMIAL5D");

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
  } //end namespace Unprec TwoFlavorRatioConvConvWilsonFermMonomialEnv


  // Constructor
  UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5D::UnprecTwoFlavorRatioConvConvWilsonTypeFermMonomial5D(
    const TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << endl;

    invParam_num = param.numer.invParam;

    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.numer.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct fermion action= " << param.numer.fermact.id << endl;

      WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
	TheWilsonTypeFermAct5DFactory::Instance().createObject(param.numer.fermact.id, fermact_reader, param.numer.fermact.path);

      UnprecWilsonTypeFermAct5D<T,P,Q>* downcast = 
	dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct5D" << endl;
	QDP_abort(1);
      }

      fermact_num = downcast;
    }

    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.denom.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct fermion action= " << param.denom.fermact.id << endl;

      WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
	TheWilsonTypeFermAct5DFactory::Instance().createObject(param.denom.fermact.id, fermact_reader, param.denom.fermact.path);

      UnprecWilsonTypeFermAct5D<T,P,Q>* downcast = 
	dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct5D" << endl;
	QDP_abort(1);
      }

      fermact_den = downcast;
    }
    //*********************************************************************

    if (fermact_num->size() != fermact_den->size()) 
    {
      QDPIO::cerr << "Error: numerator action has to have the same length in the 5th dimension as the denominator action." << endl;
      QDPIO::cerr << "N5 in FermionActionNum " << fermact_num->size() << endl;
      QDPIO::cerr << "N5 in FermionActionDen " << fermact_den->size() << endl;
      QDP_abort(1);
    }
      

    //*********************************************************************
    // Get Chronological predictor
    {
      AbsChronologicalPredictor5D<LatticeFermion>* tmp = 0x0;
      if( param.predictor.xml == "" ) {
	// No predictor specified use zero guess
	tmp = new ZeroGuess5DChronoPredictor(fermact_num->size());
      }
      else 
      {
	try 
	{ 
	  std::istringstream chrono_is(param.predictor.xml);
	  XMLReader chrono_xml(chrono_is);
	  tmp = The5DChronologicalPredictorFactory::Instance().createObject(param.predictor.id, 
									    fermact_num->size(),
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
    }
    //*********************************************************************

    QDPIO::cout << "Initing PF field" << endl;
    getPhi().resize( fermact_num->size() );

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }

} //end namespace Chroma


