// $Id: unprec_one_flavor_ratio_rat_conv_monomial_w.cc,v 3.2 2008-07-22 17:41:55 bjoo Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "update/molecdyn/monomial/unprec_one_flavor_ratio_rat_conv_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma 
{ 
 
  namespace UnprecOneFlavorWilsonTypeFermRatioRatConvMonomialEnv 
  {
    namespace
    {
      //! Does the work
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path)
      {
	return new UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial(
	  OneFlavorWilsonTypeFermRatioRatConvMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name("ONE_FLAVOR_UNPREC_FERM_RATIO_RAT_CONV_MONOMIAL");

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActs4DEnv::registerAll();
	success &= RationalApproxAggregateEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  } //end namespace Unprec OneFlavorWilsonFermRatioRatConvMonomialEnv


  // Constructor
  UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial::UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial(
    const OneFlavorWilsonTypeFermRatioRatConvMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << endl;

    num_pf             = param.num_pf;
    actionInvParam_num = param.numer.action.invParam;
    forceInvParam_num  = param.numer.force.invParam;

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

#if 0
      WilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial()" << endl;
	QDP_abort(1);
      }
#endif
      fermact_num = tmp_act;    
    }

    //*********************************************************************
    // Action rational approx
    {
      std::istringstream is(param.numer.action.ratApprox.xml);
      XMLReader approx_reader(is);
      QDPIO::cout << "Construct action rational approx= " << param.numer.action.ratApprox.id << endl;

      Handle<RationalApprox> approx(TheRationalApproxFactory::Instance().createObject(
				      param.numer.action.ratApprox.id, 
				      approx_reader, 
				      param.numer.action.ratApprox.path));

      (*approx)(spfe_num, sipfe_num);
    }

    //*********************************************************************
    // Force rational approx
    {
      std::istringstream is(param.numer.force.ratApprox.xml);
      XMLReader approx_reader(is);
      QDPIO::cout << "Construct force rational approx= " << param.numer.force.ratApprox.id << endl;

      Handle<RationalApprox> approx(TheRationalApproxFactory::Instance().createObject(
				      param.numer.force.ratApprox.id, 
				      approx_reader, 
				      param.numer.force.ratApprox.path));

      RemezCoeff_t  fipfe_num;  // discard
      (*approx)(fpfe_num, fipfe_num);
    }
    //*********************************************************************

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

#if 0
      WilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial()" << endl;
	QDP_abort(1);
      }
#endif
      fermact_den = tmp_act;    
    }
    //*********************************************************************

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }


} //end namespace Chroma


