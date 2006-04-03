// $Id: sink_smearing_aggregate.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief All make sink constructors
 */

#include "meas/sinks/sink_smearing_aggregate.h"

#include "meas/sinks/pt_sink_smearing.h"
#include "meas/sinks/sh_sink_smearing.h"

namespace Chroma
{

  //! Registration aggregator
  namespace QuarkSinkSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sinks
      success &= PointQuarkSinkSmearingEnv::registered;
      success &= ShellQuarkSinkSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
