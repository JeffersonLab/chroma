// -*- C++ -*-
// $Id: fermbcs_reader_w.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Fermionic boundary condition reader
 */

#ifndef __fermbcs_reader_w_h__
#define __fermbcs_reader_w_h__

#include "fermbc.h"
#include "handle.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup fermbc */
  namespace WilsonTypeFermBCEnv
  {
    //! Helper function for the FermionAction readers
    /*! 
     * \ingroup fermbc
     *
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC<LatticeFermion> > reader(XMLReader& xml_in, const std::string& path);
  }


  //! Registration aggregator
  /*! \ingroup fermbc */
  namespace WilsonTypeFermBCArrayEnv
  {
    //! Helper function for the FermionAction readers
    /*! 
     * \ingroup fermbc
     *
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC< multi1d<LatticeFermion> > > reader(XMLReader& xml_in, const std::string& path);
  }
}

#endif
