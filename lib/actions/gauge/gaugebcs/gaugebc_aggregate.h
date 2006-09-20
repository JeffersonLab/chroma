// -*- C++ -*-
// $Id: gaugebc_aggregate.h,v 3.1 2006-09-20 20:28:01 edwards Exp $
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
    bool registerAll();

    //! Helper function for the GaugeAction readers
    Handle<GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
											const std::string& path);
  }
}

#endif
