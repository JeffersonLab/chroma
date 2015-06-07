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
    static bool registered = false;

    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= StaggeredTypeSimpleFermBCEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
