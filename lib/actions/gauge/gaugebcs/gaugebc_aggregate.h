// -*- C++ -*-
// $Id: gaugebc_aggregate.h,v 1.1 2005-01-13 02:27:42 edwards Exp $
/*! \file
 *  \brief All gauge boundary conditions
 */

#ifndef __gaugebcs_h__
#define __gaugebcs_h__

#include "chromabase.h"
#include "handle.h"
#include "gaugebc.h"

namespace Chroma
{
  //! Registration aggregator
  namespace GaugeTypeGaugeBCEnv
  {
    extern const bool registered;

    //! Helper function for the GaugeAction readers
    Handle<GaugeBC> reader(XMLReader& xml_in, const std::string& path);
  }
}

#endif
