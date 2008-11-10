// $Id: eoprec_constdet_two_flavor_ratio_conv_conv_monomial_w.cc,v 3.2 2008-11-10 17:59:07 bjoo Exp $
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

    invParam_num = param.numer.invParam;

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
    }
    //*********************************************************************

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }

} //end namespace Chroma


