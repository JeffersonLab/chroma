// $Id: eoprec_constdet_one_flavor_ratio_rat_rat_monomial5d_w.cc,v 3.1 2008-05-23 21:31:32 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/eoprec_constdet_one_flavor_ratio_rat_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatRatMonomial5DEnv 
  {
    namespace
    {
      //! Callback
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path)
      {
	return new EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatRatMonomial5D(
	  OneFlavorWilsonTypeFermRatioRatRatMonomialParams(xml, path));
      }
    
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("ONE_FLAVOR_EOPREC_CONSTDET_FERM_RATIO_RAT_RAT_MONOMIAL5D");

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
  } //end namespace EvenOddPrecConstDet OneFlavorWilsonFermRatioRatRatMonomialEnv


  // Constructor
  EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatRatMonomial5D::EvenOddPrecConstDetOneFlavorWilsonTypeFermRatioRatRatMonomial5D(
    const OneFlavorWilsonTypeFermRatioRatRatMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << endl;

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
      QDPIO::cout << "Construct fermion action= " << param.numer.fermact.id << endl;

      WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
	TheWilsonTypeFermAct5DFactory::Instance().createObject(param.numer.fermact.id, 
							     fermact_reader, 
							     param.numer.fermact.path);

      EvenOddPrecWilsonTypeFermAct5D<T,P,Q>* downcast=dynamic_cast<EvenOddPrecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);
      
      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct5D to EvenOddPrecWilsonTypeFermAct5D" << endl;
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
	TheWilsonTypeFermAct5DFactory::Instance().createObject(param.denom.fermact.id, 
							     fermact_reader, 
							     param.denom.fermact.path);

      EvenOddPrecWilsonTypeFermAct5D<T,P,Q>* downcast=dynamic_cast<EvenOddPrecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) 
      {
	QDPIO::cerr << __func__ << ": unable to downcast FermAct5D to EvenOddPrecWilsonTypeFermAct5D" << endl;
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

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }

} //end namespace Chroma


