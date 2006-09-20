// $Id: spin_insertion_aggregate.cc,v 1.2 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief All spin insertion constructors
 */

#include "meas/hadron/spin_insertion_aggregate.h"

#include "meas/hadron/no_spin_insertion.h"
#include "meas/hadron/simple_spin_insertion_w.h"

namespace Chroma
{

  // Registration aggregator
  namespace SpinInsertionEnv
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
	success &= NoSpinInsertionEnv::registerAll();
	success &= SimpleSpinInsertionEnv::registerAll();
	registered = true;
      }
      return success;
    }
  }

}
