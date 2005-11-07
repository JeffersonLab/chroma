// $Id: source_const_aggregate.cc,v 2.1 2005-11-07 06:30:06 edwards Exp $
/*! \file
 *  \brief All make source constructors
 */

#include "meas/sources/source_const_aggregate.h"

#include "meas/sources/pt_source_const.h"
#include "meas/sources/sh_source_const.h"

namespace Chroma
{

  //! Registration aggregator
  namespace PropSourceConstructionEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PointPropSourceConstEnv::registered;
      success &= ShellPropSourceConstEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }


  //! Registration aggregator
  namespace FermSourceConstructionEnv
  {
    bool registerAll() 
    {
      bool success = true;

      // Sources
      success &= PointFermSourceConstEnv::registered;
      success &= ShellFermSourceConstEnv::registered;

      return success;
    }

    const bool registered = registerAll();
  }

}
