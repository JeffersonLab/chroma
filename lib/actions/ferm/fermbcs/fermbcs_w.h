// -*- C++ -*-
// $Id: fermbcs_w.h,v 2.0 2005-09-25 21:04:27 edwards Exp $
/*! \file
 *  \brief All fermionic boundary conditions
 */

#ifndef __fermbcs_w_h__
#define __fermbcs_w_h__

#include "chromabase.h"

namespace Chroma
{
  //! Registration aggregator
  namespace WilsonTypeFermBCEnv
  {
    extern const bool registered;

    //! Helper function for the FermionAction readers
    /*! 
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC<LatticeFermion> > reader(XMLReader& xml_in, const std::string& path);
  }


  //! Registration aggregator
  namespace WilsonTypeFermBCArrayEnv
  {
    extern const bool registered;

    //! Helper function for the FermionAction readers
    /*! 
     * This structure should not be replicated. This routine helps maintain
     * backwards compatibility with the FermionAction readers by looking for
     * either the "boundary" tag or the FermionBC group
     */
    Handle< FermBC< multi1d<LatticeFermion> > > reader(XMLReader& xml_in, const std::string& path);
  }
}

#endif
