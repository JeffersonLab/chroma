/*! \file
 *  \brief All make sink constructors
 */

#include "meas/sinks/sink_smearing_aggregate.h"

#include "meas/sinks/pt_sink_smearing.h"
#include "meas/sinks/sh_sink_smearing.h"
#include "meas/sinks/norm_sh_sink_smearing.h"
#include "meas/sinks/wall_sink_smearing.h"

namespace Chroma
{

  //! Registration aggregator
  namespace QuarkSinkSmearingEnv
  {
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Sinks
	success &= PointQuarkSinkSmearingEnv::registerAll();
	success &= ShellQuarkSinkSmearingEnv::registerAll();
	success &= NormShellQuarkSinkSmearingEnv::registerAll();
	success &= WallQuarkSinkSmearingEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
