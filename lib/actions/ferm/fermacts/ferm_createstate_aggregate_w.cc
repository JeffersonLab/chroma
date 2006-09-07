// $Id: ferm_createstate_aggregate_w.cc,v 3.4 2006-09-07 04:24:28 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#include "chromabase.h"

#include "actions/ferm/fermbcs/fermbcs_aggregate_w.h"

#include "actions/ferm/fermacts/ferm_createstate_aggregate_w.h"
#include "actions/ferm/fermacts/simple_fermstate_w.h"
#include "actions/ferm/fermacts/stout_fermstate_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace CreateFermStateEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // All ferm bcs
      success &= WilsonTypeFermBCEnv::registered;

      // All fermstates
      success &= CreateSimpleFermStateEnv::registered;
      success &= CreateStoutFermStateEnv::registered;
      success &= CreateSLICFermStateEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
