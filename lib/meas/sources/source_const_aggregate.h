// -*- C++ -*-
// $Id: source_const_aggregate.h,v 2.1 2005-11-07 06:30:06 edwards Exp $
/*! \file
 *  \brief All make source constructors
 */

#ifndef __source_const_aggregate_h__
#define __source_const_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup sources */
  namespace PropSourceConstructionEnv
  {
    extern const bool registered;
  }

  //! Registration aggregator
  /*! @ingroup sources */
  namespace FermSourceConstructionEnv
  {
    extern const bool registered;
  }
}

#endif
