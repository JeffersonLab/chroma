// -*- C++ -*-
// $Id: source_construction.h,v 1.1 2005-10-28 21:06:41 edwards Exp $

/*! @file
 * @brief Source construction
 */

#ifndef __source_construction_h__
#define __source_construction_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for source construction
  /*! @ingroup sources
   *
   * Supports creation and application of sources
   */
  template<typename T>
  class SourceConstruction
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SourceConstruction() {}

    //! Construct the source
    virtual T operator()(const multi1d<LatticeColorMatrix>& u) const = 0;
  };

}


#endif
