// $Id: monomial_aggregate_w.cc,v 2.4 2006-01-16 02:10:04 bjoo Exp $
/*! \file
 *  \brief Fermion monomial aggregator
 */

#include "update/molecdyn/monomial/monomial_aggregate_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_monomial_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_hasenbusch_monomial_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/prec_constdet_two_flavor_monomial5d_w.h"

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/prec_constdet_one_flavor_rat_monomial_w.h"

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/prec_constdet_one_flavor_rat_monomial5d_w.h"



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
    
      return success;
    }

    const bool registered = registerAll();
  }

}
