// -*- C++ -*-
/*! \file
 *  \brief All gauge field initializers
 */

#ifndef __gauge_init_aggregate_w_h__
#define __gauge_init_aggregate_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup gauge */
  namespace GaugeInitEnv
  {
    bool registerAll();
  }
}

#endif
