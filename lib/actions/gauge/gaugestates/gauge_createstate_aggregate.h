// -*- C++ -*-
/*! \file
 *  \brief All gauge create-state method
 */

#ifndef __gauge_createstate_aggregate_h__
#define __gauge_createstate_aggregate_h__

#include "chromabase.h"
#include "create_state.h"
#include "io/xml_group_reader.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup gaugestates */
  namespace CreateGaugeStateEnv
  {
    bool registerAll();

    //! Helper function for the CreateGaugeState readers
    Handle< CreateGaugeState< multi1d<LatticeColorMatrix>, 
			      multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								     const std::string& path);

    //! Returns a simple createstate group
    GroupXML_t   nullXMLGroup();
  }

}

#endif
