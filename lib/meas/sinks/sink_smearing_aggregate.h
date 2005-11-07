// -*- C++ -*-
// $Id: sink_smearing_aggregate.h,v 1.1 2005-11-07 06:24:09 edwards Exp $
/*! \file
 *  \brief All make sink constructors
 */

#ifndef __sink_smearing_aggregate_h__
#define __sink_smearing_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup sinks */
  namespace PropSinkSmearingEnv
  {
    extern const bool registered;
  }

  //! Registration aggregator
  /*! @ingroup sinks */
  namespace FermSinkSmearingEnv
  {
    extern const bool registered;
  }
}

#endif
