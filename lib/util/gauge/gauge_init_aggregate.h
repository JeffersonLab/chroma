// -*- C++ -*-
// $Id: gauge_init_aggregate.h,v 3.1 2007-02-04 22:06:42 edwards Exp $
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
