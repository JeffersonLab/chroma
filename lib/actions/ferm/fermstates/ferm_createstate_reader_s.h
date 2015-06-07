// -*- C++ -*-
/*! \file
 *  \brief All ferm create-state method
 */

#ifndef __ferm_createstate_aggregate_s_h__
#define __ferm_createstate_aggregate_s_h__

#include "chromabase.h"
#include "create_state.h"

namespace Chroma
{
  //! Registration aggregator
  /*! @ingroup fermstates */
  namespace StaggeredCreateFermStateEnv
  {
    //extern const bool registered;

    //! Helper function for the CreateFermState readers
    Handle< CreateFermState< LatticeStaggeredFermion,
			     multi1d<LatticeColorMatrix>, 
			     multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								    const std::string& path);
  }

}

#endif
