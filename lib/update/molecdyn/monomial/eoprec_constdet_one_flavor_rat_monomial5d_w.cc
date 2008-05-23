// $Id: eoprec_constdet_one_flavor_rat_monomial5d_w.cc,v 3.2 2008-05-23 21:31:32 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/eoprec_constdet_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    namespace
    {
      //! Callback
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path)
      {
	return new EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(
	  OneFlavorWilsonTypeFermRatMonomialParams(xml, path));
      }
    
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("ONE_FLAVOR_EOPREC_CONSTDET_FERM_RAT_MONOMIAL5D");

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
  } //end namespace EvenOddPrecConstDet OneFlavorWilsonFermRatMonomialEnv


  // Constructor
  EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D::EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(
    const OneFlavorWilsonTypeFermRatMonomialParams& param) 
  {
    START_CODE();

    QDPIO::cout << "Constructor: " << __func__ << endl;

    actionInvParam = param.numer.action.invParam;
    forceInvParam  = param.numer.force.invParam;
    num_pf         = param.num_pf;

    QDPIO::cout << "Number of pseudoferms= " << num_pf << endl;
    //*********************************************************************
    // Fermion action
    {
      std::istringstream is(param.numer.fermact.xml);
      XMLReader fermact_reader(is);
      QDPIO::cout << "Construct fermion action= " << param.numer.fermact.id << endl;
    
      WilsonTypeFermAct5D<T,P,Q>* tmp_act = TheWilsonTypeFermAct5DFactory::Instance().createObject(
	param.numer.fermact.id, 
	fermact_reader, 
	param.numer.fermact.path);

      EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>* downcast = 
	dynamic_cast<EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

      // Check success of the downcast 
      if( downcast == 0x0 ) {
	QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecConstDetWilsonTypeFermAct5D in EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D()" << endl;
	QDP_abort(1);
      }

      fermact = downcast;    
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

      (*approx)(spfe, sipfe);
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

      RemezCoeff_t  fipfe;  // discard
      (*approx)(fpfe, fipfe);
    }
    //*********************************************************************

    QDPIO::cout << "Finished constructing: " << __func__ << endl;
    
    END_CODE();
  }


} //end namespace Chroma


