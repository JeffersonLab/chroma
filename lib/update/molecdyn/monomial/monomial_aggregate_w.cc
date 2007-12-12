// $Id: monomial_aggregate_w.cc,v 3.7 2007-12-12 21:42:58 bjoo Exp $
/*! \file
 *  \brief Fermion monomial aggregator
 */

#include "update/molecdyn/monomial/monomial_aggregate_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/eoprec_logdet_two_flavor_monomial_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_hasenbusch_monomial_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_monomial5d_w.h"

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_one_flavor_rat_monomial_w.h"

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_one_flavor_rat_monomial5d_w.h"

//#include "update/molecdyn/monomial/unprec_two_flavor_polynomial_monomial_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_polynomial_monomial_w.h"

//#include "update/molecdyn/monomial/unprec_two_flavor_polyprec_monomial_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_polyprec_monomial_w.h"

#include "update/molecdyn/monomial/eoprec_logdet_ee_monomial_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial5d_w.h"
#include "update/molecdyn/monomial/eoprec_constdet_two_flavor_hasenbusch_monomial5d_w.h"

#include "update/molecdyn/monomial/fixed_random_ferm_monomial.h"
#include "update/molecdyn/monomial/central_tprec_logdet_tt_monomial_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeFermMonomialAggregrateEnv
  {
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// 4D Ferm Monomials
	success &= UnprecTwoFlavorWilsonTypeFermMonomialEnv::registerAll();
	success &= EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomialEnv::registerAll();
	success &= EvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomialEnv::registerAll();
   
	// 4D Ferm Monomials
	success &= UnprecOneFlavorWilsonTypeFermRatMonomialEnv::registerAll();
	success &= EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomialEnv::registerAll();
    
	// 5D Ferm Monomials
	success &= UnprecTwoFlavorWilsonTypeFermMonomial5DEnv::registerAll();
	success &= EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5DEnv::registerAll();
    
	// 5D Ferm Monomials
	success &= UnprecOneFlavorWilsonTypeFermRatMonomial5DEnv::registerAll();
	success &= EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5DEnv::registerAll();

	// Hasenbusch Monomials
	success &= UnprecTwoFlavorHasenbuschWilsonTypeFermMonomialEnv::registerAll();
	success &= EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomialEnv::registerAll();
    
	// Polynomial preconditioning Monomials
//      success &= UnprecTwoFlavorPolynomialWilsonTypeFermMonomialEnv::registerAll();
	success &= EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomialEnv::registerAll();
//      success &= UnprecTwoFlavorPolyPrecWilsonTypeFermMonomialEnv::registerAll();
	success &= EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomialEnv::registerAll();

	// Even Even part of a logdet monomial
	success &= EvenOddPrecLogDetEvenEvenMonomial4DEnv::registerAll();

	// 5D Hasenbusch Monomials
	success &= EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomial5DEnv::registerAll();

	success &= UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial5DEnv::registerAll();

	success &= FixedRandomFermMonomial4DEnv::registerAll();

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4
	success &=  CentralTimePrecLogDetTTMonomial4DEnv::registerAll();
#endif
#endif
#endif
	registered = true;
      }
      return success;
    }
  }

}
