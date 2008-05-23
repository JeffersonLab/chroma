// $Id: unprec_two_flavor_ratio_conv_rat_monomial_w.cc,v 3.1 2008-05-23 21:31:36 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_ratio_conv_rat_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorRatioConvRatWilsonTypeFermMonomialEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
      {
	return new UnprecTwoFlavorRatioConvRatWilsonTypeFermMonomial(
	  TwoFlavorRatioConvRatWilsonTypeFermMonomialParams(xml, path));
      }
 
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("TWO_FLAVOR_UNPREC_RATIO_CONV_RAT_FERM_MONOMIAL");

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
  } //end namespace Unprec TwoFlavorRatioConvRatWilsonFermMonomialEnv


  // Constructor
  UnprecTwoFlavorRatioConvRatWilsonTypeFermMonomial::UnprecTwoFlavorRatioConvRatWilsonTypeFermMonomial(
    const TwoFlavorRatioConvRatWilsonTypeFermMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << endl;

    invParam_num       = param.numer.invParam;
    actionInvParam_den = param.denom.action.invParam;
    forceInvParam_den  = param.denom.force.invParam;

    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.numer.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct fermion action= " << param.numer.fermact.id << endl;

      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.numer.fermact.id, 
							     fermact_reader, 
							     param.numer.fermact.path);

      UnprecWilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct" << endl;
	QDP_abort(1);
      }

      fermact_num = downcast;    
    }
    //*********************************************************************

    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.denom.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct fermion action= " << param.denom.fermact.id << endl;

      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.denom.fermact.id, 
							     fermact_reader, 
							     param.denom.fermact.path);

      UnprecWilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct" << endl;
	QDP_abort(1);
      }

      fermact_den = downcast;    
    }

    //*********************************************************************
    // Action rational approx
    {
      std::istringstream is(param.denom.action.ratApprox.xml);
      XMLReader approx_reader(is);
      QDPIO::cout << "Construct action rational approx= " << param.denom.action.ratApprox.id << endl;

      Handle<RationalApprox> approx(TheRationalApproxFactory::Instance().createObject(
				      param.denom.action.ratApprox.id, 
				      approx_reader, 
				      param.denom.action.ratApprox.path));

      (*approx)(spfe_den, sipfe_den);
    }

    //*********************************************************************
    // Force rational approx
    {
      std::istringstream is(param.denom.force.ratApprox.xml);
      XMLReader approx_reader(is);
      QDPIO::cout << "Construct force rational approx= " << param.denom.force.ratApprox.id << endl;

      Handle<RationalApprox> approx(TheRationalApproxFactory::Instance().createObject(
				      param.denom.force.ratApprox.id, 
				      approx_reader, 
				      param.denom.force.ratApprox.path));

      RemezCoeff_t  fipfe_den;  // discard
      (*approx)(fpfe_den, fipfe_den);
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
	QDPIO::cerr << "Failed to create ZeroGuess4DChronoPredictor" << endl;
	QDP_abort(1);
      }

      chrono_predictor = tmp;
    }
    //*********************************************************************

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }

} //end namespace Chroma


