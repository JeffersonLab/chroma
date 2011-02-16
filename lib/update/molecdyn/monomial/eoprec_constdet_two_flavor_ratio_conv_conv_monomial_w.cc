// $Id: eoprec_constdet_two_flavor_ratio_conv_conv_monomial_w.cc,v 3.3 2009-06-04 20:30:27 bjoo Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_ratio_conv_conv_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"

#include "update/molecdyn/predictor/null_predictor.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomialEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
      {
	return new EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial(
	  TwoFlavorRatioConvConvWilsonTypeFermMonomialParams(xml, path));
      }
    
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("TWO_FLAVOR_EOPREC_CONSTDET_RATIO_CONV_CONV_FERM_MONOMIAL");

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
  EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial::EvenOddPrecConstDetTwoFlavorRatioConvConvWilsonTypeFermMonomial(
    const TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << endl;

    if ( param.numer.invParam.id == "NULL" ) { 
      QDPIO::cerr << "Numerator inverter parameter is NULL" << endl;
      QDP_abort(1);
    }
    invParam_num = param.numer.invParam;

    if( param.denom.invParam.id == "NULL" ) { 
      QDPIO::cerr << "WARNING: Denominator inverter parameter is NULL." << endl;
      QDPIO::cerr << "WARNING: Assuming it is same as numerator " << endl;
      invParam_den = param.numer.invParam;
    }
    else {
      invParam_den = param.denom.invParam;
    }

    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.numer.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct numer fermion action= " << param.numer.fermact.id << endl;

      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.numer.fermact.id, 
							     fermact_reader, 
							     param.numer.fermact.path);

      EvenOddPrecWilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<EvenOddPrecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to EvenOddPrecConstDetWilsonTypeFermAct" << endl;
	QDP_abort(1);
      }

      fermact_num = downcast;    
    }

    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.denom.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct denom fermion action= " << param.denom.fermact.id << endl;

      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.denom.fermact.id, 
							     fermact_reader, 
							     param.denom.fermact.path);

      EvenOddPrecWilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<EvenOddPrecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to EvenOddPrecConstDetWilsonTypeFermAct" << endl;
	QDP_abort(1);
      }

      fermact_den = downcast;    
    }
    //*********************************************************************

    //*********************************************************************
    // Get Chronological predictor
    {
      AbsChronologicalPredictor4D<LatticeFermion>* tmp = 0x0;
      if( param.predictor.xml == "" ) {
	// No predictor specified use: Null guess
	QDPIO::cout << "No predictor specified. Using NULL Predictor" << endl;
	tmp = new Null4DChronoPredictor();
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
    }
    //*********************************************************************

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }

} //end namespace Chroma


