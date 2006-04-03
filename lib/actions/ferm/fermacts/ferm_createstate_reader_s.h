// -*- C++ -*-
// $Id: ferm_createstate_reader_s.h,v 3.0 2006-04-03 04:58:45 edwards Exp $
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
  /*! @ingroup fermacts */
  namespace StaggeredCreateFermStateEnv
  {
    extern const bool registered;

    //! Helper function for the CreateFermState readers
    Handle< CreateFermState< LatticeStaggeredFermion,
			     multi1d<LatticeColorMatrix>, 
			     multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								    const std::string& path);
  }

}

#endif
