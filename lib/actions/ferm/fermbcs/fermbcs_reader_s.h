// -*- C++ -*-
// $Id: fermbcs_reader_s.h,v 3.1 2006-08-18 15:52:43 edwards Exp $
/*! \file
 *  \brief Fermionic boundary condition reader
 */

#ifndef __fermbcs_reader_s_h__
#define __fermbcs_reader_s_h__

#include "fermbc.h"
#include "handle.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup fermbcs */
  namespace StaggeredTypeFermBCEnv
  {
    //! Helper function for the FermionAction readers
    /*!
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC<LatticeStaggeredFermion,
                   multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
							  const std::string& path);
  }
}

#endif
