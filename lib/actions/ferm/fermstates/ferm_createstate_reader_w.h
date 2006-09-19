// -*- C++ -*-
// $Id: ferm_createstate_reader_w.h,v 1.1 2006-09-19 17:53:37 edwards Exp $
/*! \file
 *  \brief All ferm create-state method
 */

#ifndef __ferm_createstate_reader_w_h__
#define __ferm_createstate_reader_w_h__

#include "create_state.h"
#include "handle.h"

namespace Chroma
{
  //! State reader
  /*! @ingroup fermstates */
  namespace CreateFermStateEnv
  {
    //! Helper function for the CreateFermState readers
    Handle< CreateFermState< LatticeFermion,
			     multi1d<LatticeColorMatrix>, 
			     multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
								    const std::string& path);
  }

}

#endif
