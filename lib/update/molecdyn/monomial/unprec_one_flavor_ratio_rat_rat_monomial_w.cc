/*! @file
 * @brief One-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "update/molecdyn/monomial/unprec_one_flavor_ratio_rat_rat_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma 
{ 
 
  namespace UnprecOneFlavorWilsonTypeFermRatioRatRatMonomialEnv 
  {
    namespace
    {
      //! Does the work
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const std::string& path)
      {
	return new UnprecOneFlavorWilsonTypeFermRatioRatRatMonomial(
	  OneFlavorWilsonTypeFermRatioRatRatMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name("ONE_FLAVOR_UNPREC_FERM_RATIO_RAT_RAT_MONOMIAL");

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
  } //end namespace Unprec OneFlavorWilsonFermRatioRatRatMonomialEnv


  // Constructor
  UnprecOneFlavorWilsonTypeFermRatioRatRatMonomial::UnprecOneFlavorWilsonTypeFermRatioRatRatMonomial(
    const OneFlavorWilsonTypeFermRatioRatRatMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << std::endl;

    num_pf             = param.num_pf;
    actionInvParam_num = param.numer.action.invParam;
    forceInvParam_num  = param.numer.force.invParam;
    actionInvParam_den = param.denom.action.invParam;
    forceInvParam_den  = param.denom.force.invParam;

    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.numer.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct fermion action= " << param.numer.fermact.id << std::endl;

      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.numer.fermact.id, 
							     fermact_reader, 
							     param.numer.fermact.path);
#if 0
      WilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct" << std::endl;
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
      QDPIO::cout << "Construct action rational approx= " << param.numer.action.ratApprox.id << std::endl;

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
      QDPIO::cout << "Construct force rational approx= " << param.numer.force.ratApprox.id << std::endl;

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
      QDPIO::cout << "Construct fermion action= " << param.denom.fermact.id << std::endl;

      WilsonTypeFermAct<T,P,Q>* tmp_act = 
	TheWilsonTypeFermActFactory::Instance().createObject(param.denom.fermact.id, 
							     fermact_reader, 
							     param.denom.fermact.path);

#if 0
      WilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct to UnprecWilsonTypeFermAct" << std::endl;
	QDP_abort(1);
      }
#endif
      fermact_den = tmp_act;    
    }

    //*********************************************************************
    // Action rational approx
    {
      std::istringstream is(param.denom.action.ratApprox.xml);
      XMLReader approx_reader(is);
      QDPIO::cout << "Construct action rational approx= " << param.denom.action.ratApprox.id << std::endl;

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
      QDPIO::cout << "Construct force rational approx= " << param.denom.force.ratApprox.id << std::endl;

      Handle<RationalApprox> approx(TheRationalApproxFactory::Instance().createObject(
				      param.denom.force.ratApprox.id, 
				      approx_reader, 
				      param.denom.force.ratApprox.path));

      RemezCoeff_t  fipfe_den;  // discard
      (*approx)(fpfe_den, fipfe_den);
    }
    //*********************************************************************

    QDPIO::cout << "Finished constructing: " << __func__ << std::endl;
    
    END_CODE();
  }


} //end namespace Chroma


