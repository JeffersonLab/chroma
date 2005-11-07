// -*- C++ -*-
// $Id: source_construction.h,v 1.2 2005-11-07 06:28:35 edwards Exp $

/*! @file
 * @brief Source construction
 */

#ifndef __source_construction_h__
#define __source_construction_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for quark source construction
  /*! @ingroup sources
   *
   * Supports creation of quark sources
   */
  template<typename T>
  class QuarkSourceConstruction
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~QuarkSourceConstruction() {}

    //! Construct the source
    virtual T operator()(const multi1d<LatticeColorMatrix>& u) const = 0;
  };

}


#endif
