// $Id: unprec_one_flavor_rat_monomial5d_w.cc,v 3.3 2006-09-20 20:28:05 edwards Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "update/molecdyn/monomial/genapprox.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma 
{ 
 
  namespace UnprecOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    namespace
    {
      //! Callback
      Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path)
      {
	QDPIO::cout << "Create Monomial: " << name << endl;

	return new UnprecOneFlavorWilsonTypeFermRatMonomial5D(
	  OneFlavorWilsonTypeFermRatMonomial5DParams(xml, path));
      }
    
      //! Local registration flag
      bool registered = false;
    }

    const std::string name("ONE_FLAVOR_UNPREC_FERM_RAT_MONOMIAL5D");

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
  } //end namespace Unprec OneFlavorWilsonFermRatMonomialEnv


  // Constructor
  UnprecOneFlavorWilsonTypeFermRatMonomial5D::UnprecOneFlavorWilsonTypeFermRatMonomial5D(
    const OneFlavorWilsonTypeFermRatMonomial5DParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;
    nthRoot   = param.nthRoot;
    nthRootPV = param.nthRoot;

    std::istringstream is(param.fermact.xml);
    XMLReader fermact_reader(is);
    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial5D: construct " << param.fermact.id << endl;

    WilsonTypeFermAct5D<T,P,Q>* tmp_act = 
      TheWilsonTypeFermAct5DFactory::Instance().createObject(param.fermact.id, fermact_reader, param.fermact.path);

    UnprecWilsonTypeFermAct5D<T,P,Q>* downcast = 
      dynamic_cast<UnprecWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct5D in UnprecOneFlavorWilsonTypeFermRatMonomial5D()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    //*********************************************************************
    // Remez approx
    // M term
    QDPIO::cout << "Normal operator PFE" << endl;
    generateApprox(fpfe, spfe, sipfe,
		   param.remez.lowerMin, param.remez.upperMax, 
		   -param.expNumPower, 2*param.expDenPower*nthRoot, 
		   param.remez.degree, param.remez.degree,
		   param.remez.digitPrecision);

    // PV term
    QDPIO::cout << "PV operator PFE" << endl;
    generateApprox(fpvpfe, spvpfe, sipvpfe,
		   param.remez.lowerMinPV, param.remez.upperMaxPV, 
		   param.expNumPower, 2*param.expDenPower*nthRootPV, 
		   param.remez.degreePV, param.remez.degreePV,
		   param.remez.digitPrecision);
    //*********************************************************************

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial5D: finished " << param.fermact.id << endl;
    
    END_CODE();
  }

} //end namespace Chroma
