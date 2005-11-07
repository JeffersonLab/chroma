// -*- C++ -*-
// $Id: source_smearing_aggregate.h,v 2.1 2005-11-07 06:30:06 edwards Exp $
/*! \file
 *  \brief All source smearing
 */

#ifndef __source_smearing_aggregate_h__
#define __source_smearing_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup sources */
  namespace PropSourceSmearingEnv
  {
    extern const bool registered;
  }

  //! Registration aggregator
  /*! @ingroup sources */
  namespace FermSourceSmearingEnv
  {
    extern const bool registered;
  }
}

#endif
