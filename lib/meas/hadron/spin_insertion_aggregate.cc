// $Id: spin_insertion_aggregate.cc,v 1.1 2006-05-24 21:09:41 edwards Exp $
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
    bool registerAll() 
    {
      bool success = true;

      success &= NoSpinInsertionEnv::registered;
      success &= SimpleSpinInsertionEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
