// -*- C++ -*-
// $Id: gauge_createstate_aggregate.h,v 3.0 2006-04-03 04:58:53 edwards Exp $
/*! \file
 *  \brief All gauge create-state method
 */

#ifndef __gauge_createstate_aggregate_h__
#define __gauge_createstate_aggregate_h__

#include "chromabase.h"
#include "create_state.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup gaugeacts */
  namespace CreateGaugeStateEnv
  {
    extern const bool registered;

    //! Helper function for the CreateGaugeState readers
    Handle< CreateGaugeState< multi1d<LatticeColorMatrix>, 
			      multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								     const std::string& path);
  }

}

#endif
