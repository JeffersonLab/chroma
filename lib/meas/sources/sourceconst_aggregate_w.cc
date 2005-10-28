// $Id: sourceconst_aggregate_w.cc,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief All make source constructors
 */

#include "meas/sources/sourceconst_aggregate_w.h"

#include "meas/sources/pt_sourceconst_w.h"
#include "meas/sources/sh_sourceconst_w.h"

namespace Chroma
{

  //! Registration aggregator
  namespace SourceConstructionEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PropPointSourceConstEnv::registered;
      success &= PropShellSourceConstEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
