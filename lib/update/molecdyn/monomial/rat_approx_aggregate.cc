/*! \file
 *  \brief Rational approximation aggregator
 */

#include "update/molecdyn/monomial/rat_approx_aggregate.h"
#include "update/molecdyn/monomial/remez_rat_approx.h"
#include "update/molecdyn/monomial/read_rat_approx.h"

namespace Chroma
{

  //! Name and registration
  namespace RationalApproxAggregateEnv
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
	success &= RemezRatApproxEnv::registerAll();
	success &= ReadRatApproxEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
