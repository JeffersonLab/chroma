// $Id: fermbcs_aggregate_s.cc,v 3.1 2006-08-18 15:52:43 edwards Exp $
/*! \file
 *  \brief All Wilson-type fermion boundary conditions
 */

#include "actions/ferm/fermbcs/fermbcs_aggregate_s.h"
#include "actions/ferm/fermbcs/simple_fermbc_s.h"

namespace Chroma
{

  //! Registration aggregator
  namespace StaggeredTypeFermBCEnv
  {
    bool registerAll(void) 
    {
      bool success; 
      success = StaggeredTypeSimpleFermBCEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
