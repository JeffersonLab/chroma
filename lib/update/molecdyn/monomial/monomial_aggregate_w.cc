// $Id: monomial_aggregate_w.cc,v 3.0 2006-04-03 04:59:08 edwards Exp $
/*! \file
 *  \brief Fermion monomial aggregator
 */

#include "update/molecdyn/monomial/monomial_aggregate_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/prec_logdet_two_flavor_monomial_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_hasenbusch_monomial_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_monomial5d_w.h"

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_one_flavor_rat_monomial_w.h"

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/prec_constdet_one_flavor_rat_monomial5d_w.h"

//#include "update/molecdyn/monomial/unprec_two_flavor_polynomial_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_polynomial_monomial_w.h"

//#include "update/molecdyn/monomial/unprec_two_flavor_polyprec_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_polyprec_monomial_w.h"

#include "update/molecdyn/monomial/prec_logdet_ee_monomial_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeFermMonomialAggregrateEnv
  {
    bool registerAll() 
    {
      bool success = true; 

      // 4D Ferm Monomials
      success &= UnprecTwoFlavorWilsonTypeFermMonomialEnv::registered;
      success &= EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomialEnv::registered;
      success &= EvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomialEnv::registered;
   
      // 4D Ferm Monomials
      success &= UnprecOneFlavorWilsonTypeFermRatMonomialEnv::registered;
      success &= EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomialEnv::registered;
    
      // 5D Ferm Monomials
      success &= UnprecTwoFlavorWilsonTypeFermMonomial5DEnv::registered;
      success &= EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial5DEnv::registered;
    
      // 5D Ferm Monomials
      success &= UnprecOneFlavorWilsonTypeFermRatMonomial5DEnv::registered;
      success &= EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5DEnv::registered;

      // Hasenbusch Monomials
      success &=   UnprecTwoFlavorHasenbuschWilsonTypeFermMonomialEnv::registered;
      success &=   EvenOddPrecConstDetTwoFlavorHasenbuschWilsonTypeFermMonomialEnv::registered;
    
      // Polynomial preconditioning Monomials
//      success &=   UnprecTwoFlavorPolynomialWilsonTypeFermMonomialEnv::registered;
      success &=   EvenOddPrecConstDetTwoFlavorPolynomialWilsonTypeFermMonomialEnv::registered;
//      success &=   UnprecTwoFlavorPolyPrecWilsonTypeFermMonomialEnv::registered;
      success &=   EvenOddPrecConstDetTwoFlavorPolyPrecWilsonTypeFermMonomialEnv::registered;

      // Even Even part of a logdet monomial
      success &=  PrecLogDetEvenEvenMonomial4DEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
