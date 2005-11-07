// $Id: source_smearing_aggregate.cc,v 2.1 2005-11-07 06:30:06 edwards Exp $
/*! \file
 *  \brief All source smearing
 */

#include "meas/sources/source_smearing_aggregate.h"

#include "meas/sources/pt_source_smearing.h"
#include "meas/sources/sh_source_smearing.h"

namespace Chroma
{

  //! Registration aggregator
  namespace PropSourceSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PointPropSourceSmearingEnv::registered;
      success &= ShellPropSourceSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }


  //! Registration aggregator
  namespace FermSourceSmearingEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PointFermSourceSmearingEnv::registered;
      success &= ShellFermSourceSmearingEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
