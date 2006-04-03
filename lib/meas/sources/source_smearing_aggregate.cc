// $Id: source_smearing_aggregate.cc,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief All source smearing
 */

#include "meas/sources/source_smearing_aggregate.h"

#include "meas/sources/pt_source_smearing.h"
#include "meas/sources/sh_source_smearing.h"

namespace Chroma
{

  //! Registration aggregator
  namespace QuarkSourceSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PointQuarkSourceSmearingEnv::registered;
      success &= ShellQuarkSourceSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
