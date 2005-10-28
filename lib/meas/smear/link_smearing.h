// -*- C++ -*-
// $Id: link_smearing.h,v 1.1 2005-10-28 21:31:04 edwards Exp $
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
    virtual multi1d<LatticeColorMatrix> operator()(const multi1d<LatticeColorMatrix>& u) const = 0;
  };

}


#endif
