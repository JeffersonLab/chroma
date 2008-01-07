// -*- C++ -*-
// $Id: dilution_scheme_aggregate.h,v 1.1 2008-01-07 15:21:26 jbulava Exp $
/*! \file
 *  \brief All dilution scheme factories
 */

#ifndef __dilution_scheme_aggregate_h__
#define __dilution_scheme_aggregate_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup hadron */
  namespace DilutionSchemeEnv
  {
    bool registerAll();
  }
}

#endif
