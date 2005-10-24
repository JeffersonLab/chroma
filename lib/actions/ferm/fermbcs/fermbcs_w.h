// -*- C++ -*-
// $Id: fermbcs_w.h,v 2.1 2005-10-24 05:51:40 edwards Exp $
/*! \file
 *  \brief All fermionic boundary conditions
 */

#ifndef __fermbcs_w_h__
#define __fermbcs_w_h__

#include "fermbc.h"

namespace Chroma
{
  //! Registration aggregator
  /*! \ingroup fermbc */
  namespace WilsonTypeFermBCEnv
  {
    /*! \ingroup fermbc */
    extern const bool registered;

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
    extern const bool registered;

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
