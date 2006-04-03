// $Id: fermbcs_aggregate_s.cc,v 3.0 2006-04-03 04:58:48 edwards Exp $
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
