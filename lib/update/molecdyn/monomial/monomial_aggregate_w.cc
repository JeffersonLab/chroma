// $Id: monomial_aggregate_w.cc,v 1.1 2005-02-07 04:11:55 edwards Exp $
/*! \file
 *  \brief Fermion monomial aggregator
 */

#include "update/molecdyn/monomial/monomial_aggregate_w.h"

#include "update/molecdyn/monomial/unprec_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/prec_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/unprec_two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/prec_two_flavor_monomial5d_w.h"

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/prec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/prec_one_flavor_rat_monomial5d_w.h"

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
      success &= EvenOddPrecTwoFlavorWilsonTypeFermMonomialEnv::registered;
    
      // 4D Ferm Monomials
      success &= UnprecOneFlavorWilsonTypeFermRatMonomialEnv::registered;
      success &= EvenOddPrecOneFlavorWilsonTypeFermRatMonomialEnv::registered;
    
      // 5D Ferm Monomials
      success &= UnprecTwoFlavorWilsonTypeFermMonomial5DEnv::registered;
      success &= EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DEnv::registered;
    
      // 5D Ferm Monomials
      success &= UnprecOneFlavorWilsonTypeFermRatMonomial5DEnv::registered;
      success &= EvenOddPrecOneFlavorWilsonTypeFermRatMonomial5DEnv::registered;
    
      return success;
    }

    const bool registered = registerAll();
  }

}
