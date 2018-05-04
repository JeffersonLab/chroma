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
    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Sources
	success &= PointQuarkSourceSmearingEnv::registerAll();
	success &= ShellQuarkSourceSmearingEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
