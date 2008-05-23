// $Id: unprec_one_flavor_ratio_rat_conv_monomial5d_w.cc,v 3.1 2008-05-23 21:31:34 edwards Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/unprec_one_flavor_ratio_rat_conv_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma 
{ 
 
  namespace UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial5DEnv 
  {
    namespace
    {
      //! Callback
      Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path)
      {
	return new UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial5D(
	  OneFlavorWilsonTypeFermRatioRatConvMonomialParams(xml, path));
      }
    
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("ONE_FLAVOR_UNPREC_FERM_RATIO_RAT_CONV_MONOMIAL5D");

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActs5DEnv::registerAll();
	success &= RationalApproxAggregateEnv::registerAll();
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  } //end namespace Unprec OneFlavorWilsonFermRatioRatConvMonomialEnv


  // Constructor
  UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial5D::UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial5D(
    const OneFlavorWilsonTypeFermRatioRatConvMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << endl;

    actionInvParam_num = param.numer.action.invParam;
    forceInvParam_num  = param.numer.force.invParam;
    num_pf             = param.num_pf;

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
	QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct5D in UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial5D()" << endl;
	QDP_abort(1);
      }

      fermact_num = downcast;    
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
      QDPIO::cout << "Construct fermion action= " << param.denom.fermact.id << endl;

      WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
	TheWilsonTypeFermAct5DFactory::Instance().createObject(param.denom.fermact.id, fermact_reader, param.denom.fermact.path);

      UnprecWilsonTypeFermAct5D<T,P,Q>* downcast = 
	dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) {
	QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct5D in UnprecOneFlavorWilsonTypeFermRatioRatConvMonomial5D()" << endl;
	QDP_abort(1);
      }

      fermact_den = downcast;    
    }
    //*********************************************************************

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }

} //end namespace Chroma
