// $Id: source_const_aggregate.cc,v 2.2 2005-11-07 22:46:34 edwards Exp $
/*! \file
 *  \brief All make source constructors
 */

#include "meas/sources/source_const_aggregate.h"

#include "meas/sources/pt_source_const.h"
#include "meas/sources/sh_source_const.h"

namespace Chroma
{

  //! Registration aggregator
  namespace QuarkSourceConstructionEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PointQuarkSourceConstEnv::registered;
      success &= ShellQuarkSourceConstEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
