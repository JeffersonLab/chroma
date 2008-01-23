// $Id: unprec_one_flavor_rat_monomial_w.cc,v 3.4 2008-01-23 18:23:36 bjoo Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "update/molecdyn/monomial/genapprox.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma 
{ 
 
  namespace UnprecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    namespace
    {
      //! Does the work
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path)
      {
	QDPIO::cout << "Create Fractional Monomial: " << name << endl;

	return new UnprecOneFlavorWilsonTypeFermRatMonomial(
	  OneFlavorWilsonTypeFermRatMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name("ONE_FLAVOR_UNPREC_FERM_RAT_MONOMIAL");

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
  } //end namespace Unprec OneFlavorWilsonFermRatMonomialEnv


  // Constructor
  UnprecOneFlavorWilsonTypeFermRatMonomial::UnprecOneFlavorWilsonTypeFermRatMonomial(
    const OneFlavorWilsonTypeFermRatMonomialParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;
    nthRoot   = param.nthRoot;

    std::istringstream is(param.fermact.xml);
    XMLReader fermact_reader(is);

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial: construct " << param.fermact.id << endl;

    WilsonTypeFermAct<T,P,Q>* tmp_act = 
      TheWilsonTypeFermActFactory::Instance().createObject(param.fermact.id, 
							   fermact_reader, 
							   param.fermact.path);

    // WilsonTypeFermAct<T,P,Q>* downcast=dynamic_cast<UnprecWilsonTypeFermAct<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    // if( downcast == 0x0 ) {
    //  QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecOneFlavorWilsonTypeFermRatMonomial()" << endl;
    //  QDP_abort(1);
    //}

    //fermact = downcast;    
    fermact = tmp_act;

    //*********************************************************************
    // Remez approx
    // M term
    QDPIO::cout << "Normal operator PFE" << endl;
    generateApprox(fpfe, spfe, sipfe,
		   param.remez.lowerMin, param.remez.upperMax, 
		   -param.expNumPower, 2*param.expDenPower*nthRoot, 
		   param.remez.forceDegree, param.remez.actionDegree,
		   param.remez.digitPrecision);
    //*********************************************************************

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial: finished " << param.fermact.id << endl;
    
    END_CODE();
  }


} //end namespace Chroma


