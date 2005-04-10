// -*- C++ -*-
// $Id: gaugebc_aggregate.h,v 1.3 2005-04-10 22:25:08 edwards Exp $
/*! \file
 *  \brief Gauge boundary condition aggregator
 */

#ifndef __gaugebc_aggregate_h__
#define __gaugebc_aggregate_h__

#include "chromabase.h"
#include "handle.h"
#include "gaugebc.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup gaugebcs */
  namespace GaugeTypeGaugeBCEnv
  {
    extern const bool registered;

    //! Helper function for the GaugeAction readers
    Handle<GaugeBC> reader(XMLReader& xml_in, const std::string& path);
  }
}

#endif
