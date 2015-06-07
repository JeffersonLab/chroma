// -*- C++ -*-
/*! @file
 * @brief Link smearing
 */

#ifndef __link_smearing_h__
#define __link_smearing_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for link smearing
  /*! @ingroup sources
   *
   * Supports creation and application of link smear operators
   */
  class LinkSmearing
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~LinkSmearing() {}

    //! Do the appropriate link smearing
    virtual void operator()(multi1d<LatticeColorMatrix>& u) const = 0;
  };

}


#endif
