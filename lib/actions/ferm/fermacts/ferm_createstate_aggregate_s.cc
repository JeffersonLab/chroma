// $Id: ferm_createstate_aggregate_s.cc,v 3.0 2006-04-03 04:58:44 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/ferm_createstate_aggregate_s.h"
#include "actions/ferm/fermacts/simple_fermstate_s.h"

namespace Chroma
{

  //! Registration aggregator
  namespace StaggeredCreateFermStateEnv
  {
    bool registerAll() 
    {
      bool success = true;

      success &= StaggeredCreateSimpleFermStateEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
