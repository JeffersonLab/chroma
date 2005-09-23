// $Id: inline_io_aggregate.cc,v 1.1 2005-09-23 03:43:09 edwards Exp $
/*! \file
 *  \brief Inline IO aggregator
 */

#include "meas/inline/io/inline_io_aggregate.h"
#include "meas/inline/io/inline_write_obj.h"
#include "meas/inline/io/inline_write_erase_obj.h"
#include "meas/inline/io/inline_read_obj.h"
#include "meas/inline/io/inline_erase_obj.h"

namespace Chroma
{

  //! Name and registration
  namespace InlineIOAggregateEnv
  {
    bool registerAll() 
    {
      bool success = true; 

      // Tasks
      success &= InlineReadNamedObjEnv::registered;
      success &= InlineWriteNamedObjEnv::registered;
      success &= InlineWriteEraseNamedObjEnv::registered;
      success &= InlineEraseNamedObjEnv::registered;
      return success;
    }

    const bool registered = registerAll();
  }

}
