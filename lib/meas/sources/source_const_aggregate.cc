// $Id: source_const_aggregate.cc,v 2.3 2005-11-08 05:29:02 edwards Exp $
/*! \file
 *  \brief All make source constructors
 */

#include "meas/sources/source_const_aggregate.h"

#include "meas/sources/pt_source_const.h"
#include "meas/sources/sh_source_const.h"
#include "meas/sources/rndz2wall_source_const.h"
#include "meas/sources/wall_source_const.h"

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
      success &= RandZ2WallQuarkSourceConstEnv::registered;
      success &= WallQuarkSourceConstEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
