/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_ratio_conv_conv_multihasen_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 

	namespace UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialEnv 
	{
		namespace
		{
			//! Callback function for the factory
			Monomial< multi1d<LatticeColorMatrix>,
				multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const std::string& path) 
				{
					return new UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(
							TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams(xml, path));
				}

			//! Local registration flag
			bool registered = false;
		}

		const std::string name("TWO_FLAVOR_UNPREC_RATIO_CONV_CONV_MULTIHASEN_FERM_MONOMIAL");

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
	} //end namespace Unprec TwoFlavorRatioConvConvWilsonFermMonomialEnv


	// Constructor
	UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial::UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial(
			const TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& param) 
	{
		START_CODE();

		QDPIO::cout << "Constructor: " << __func__ << std::endl;

		if( param.numer.invParam.id == "NULL" ) {
			QDPIO::cerr << "Numerator inverter params are NULL" << std::endl;
			QDP_abort(1);
		}

		invParam_num = param.numer.invParam;

		if( param.denom.invParam.id == "NULL" ) {
			QDPIO::cerr << "WARNING: No inverter params provided for denominator." << std::endl;
			QDPIO::cerr << "WARNING: Assuming same as for numerator " << std::endl;
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
			QDPIO::cout << "Construct numer fermion action= " << param.numer.fermact.id << std::endl;

			WilsonTypeFermAct<T,P,Q>* tmp_act = 
				TheWilsonTypeFermActFactory::Instance().createObject(param.numer.fermact.id, 
						fermact_reader, 
						param.numer.fermact.path);

			UnprecWilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

			// Check success of the downcast 
			if( downcast == 0x0 ) 
			{
				QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct" << std::endl;
				QDP_abort(1);
			}

			fermact_num = downcast;    
		}

		//*********************************************************************
		// Fermion action
		{
			std::istringstream is(param.denom.fermact.xml);
			XMLReader fermact_reader(is);
			QDPIO::cout << "Construct denom fermion action= " << param.denom.fermact.id << std::endl;

			WilsonTypeFermAct<T,P,Q>* tmp_act = 
				TheWilsonTypeFermActFactory::Instance().createObject(param.denom.fermact.id, 
						fermact_reader, 
						param.denom.fermact.path);

			UnprecWilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

			// Check success of the downcast 
			if( downcast == 0x0 ) 
			{
				QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct" << std::endl;
				QDP_abort(1);
			}

			fermact_den = downcast;    
		}
		//*********************************************************************

		//*********************************************************************
		// Number of Hasenbusch term
		numHasenTerms = param.masses.size();
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
					QDPIO::cerr << "Caught Exception Reading XML: " << e << std::endl;
					QDP_abort(1);
				}
			}

			if( tmp == 0x0 ) { 
				QDPIO::cerr << "Failed to create the 4D ChronoPredictor" << std::endl;
				QDP_abort(1);
			}
			chrono_predictor = tmp;
		}
		//*********************************************************************

		QDPIO::cout << "Finished constructing: " << __func__ << std::endl;

		END_CODE();
	}
	// Sum over all Hasenbusch term to get total action
	Double UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial:: S(const AbsFieldState<P,Q>& s){
		Double S=0;
		for(int i=0; i<numHasenTerms; ++i){
			S += multihasenMonomial[i]->S(s);
		}
		return S;
	}

	// Sum over all Force terms
	void UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial::
		dsdq(P& F, const AbsFieldState<P,Q>& s){
			P F_t = zero;
			for(int i=0; i<numHasenTerms; ++i){
				numltihasenMonomial[i]->dsdq(F, s);
				F_t += F;
			}
			F = F_t;
		}

	// Refresh pseudofermion field of all Hasenbusch term
	void UnprecTwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomial::
		refreshInteralFields(const AbsFieldState<P,Q>& field_state){
			for(int i=0; i<numHasenTerms; ++i){
				multihasenMonomial[i]->refreshIntermalFields(field_state)



} //end namespace Chroma
