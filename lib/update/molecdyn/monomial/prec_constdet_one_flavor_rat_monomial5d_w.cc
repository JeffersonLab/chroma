// $Id: prec_constdet_one_flavor_rat_monomial5d_w.cc,v 3.3 2006-09-20 20:28:05 edwards Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/prec_constdet_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "update/molecdyn/monomial/genapprox.h"

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
	QDPIO::cout << "Create Monomial: " << name << endl;
	return new EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(
	  OneFlavorWilsonTypeFermRatMonomial5DParams(xml, path));
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
    const OneFlavorWilsonTypeFermRatMonomial5DParams& param) 
  {
    START_CODE();

    inv_param = param.inv_param;
    nthRoot   = param.nthRoot;
    nthRootPV = param.nthRootPV;

    std::istringstream is(param.fermact.xml);
    XMLReader fermact_reader(is);
    
    WilsonTypeFermAct5D<T,P,Q>* tmp_act = TheWilsonTypeFermAct5DFactory::Instance().createObject(
      param.fermact.id, 
      fermact_reader, 
      param.fermact.path);

    EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>* downcast = 
      dynamic_cast<EvenOddPrecConstDetWilsonTypeFermAct5D<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecConstDetWilsonTypeFermAct5D in EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D()" << endl;
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

    QDPIO::cout << "EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D: " << param.fermact.id << endl;

    END_CODE();
  }


} //end namespace Chroma


