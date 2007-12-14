// -*- C++ -*-
// $Id: dilution_operator_aggregate.h,v 1.1 2007-12-14 06:53:42 edwards Exp $
/*! \file
 *  \brief All dilution operator factories
 */

#ifndef __dilution_operator_aggregate_h__
#define __dilution_operator_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup hadron */
  namespace DilutionOperatorEnv
  {
    bool registerAll();
  }
}

#endif
