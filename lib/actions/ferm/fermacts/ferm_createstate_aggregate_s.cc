// $Id: ferm_createstate_aggregate_s.cc,v 3.1 2006-08-18 15:51:55 edwards Exp $
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
