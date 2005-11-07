// $Id: sink_smearing_aggregate.cc,v 1.1 2005-11-07 06:24:09 edwards Exp $
/*! \file
 *  \brief All make sink constructors
 */

#include "meas/sinks/sink_smearing_aggregate.h"

#include "meas/sinks/pt_sink_smearing.h"
#include "meas/sinks/sh_sink_smearing.h"

namespace Chroma
{

  //! Registration aggregator
  namespace PropSinkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sinks
      success &= PointPropSinkSmearingEnv::registered;
      success &= ShellPropSinkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

  //! Registration aggregator
  namespace FermSinkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sinks
      success &= PointFermSinkSmearingEnv::registered;
      success &= ShellFermSinkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
