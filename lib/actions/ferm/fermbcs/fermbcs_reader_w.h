// -*- C++ -*-
// $Id: fermbcs_reader_w.h,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief Fermionic boundary condition reader
 */

#ifndef __fermbcs_reader_w_h__
#define __fermbcs_reader_w_h__

#include "fermbc.h"
#include "handle.h"

namespace Chroma
{
  //! FermBC reader
  /*! \ingroup fermbcs */
  namespace WilsonTypeFermBCEnv
  {
    //! Helper function for the FermionAction readers
    /*! 
     * \ingroup fermbcs
     *
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC<LatticeFermion,
		   multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> > > reader(XMLReader& xml_in, 
							  const std::string& path);
  }

}

#endif
