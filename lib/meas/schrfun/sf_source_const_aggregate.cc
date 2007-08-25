// $Id: sf_source_const_aggregate.cc,v 1.1 2007-08-25 04:07:41 edwards Exp $
/*! \file
 *  \brief All source constructors for Schroedinger functional
 */

#include "meas/schrfun/sf_source_const_aggregate.h"

#include "meas/schrfun/sf_pt_source_const.h"
#include "meas/schrfun/sf_sh_source_const.h"
#include "meas/schrfun/sf_wall_source_const.h"

namespace Chroma
{

  //! Registration aggregator
  namespace SFSourceConstructionEnv
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
	success &= PointSFSourceConstEnv::registerAll();
	success &= ShellSFSourceConstEnv::registerAll();
	success &= WallSFSourceConstEnv::registerAll();

	registered = true;
      }
      return success;
    }
  }

}
