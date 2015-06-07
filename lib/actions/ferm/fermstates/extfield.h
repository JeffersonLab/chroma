// -*- C++ -*-

/*! @file
 * @brief External field
 */

#ifndef __extfield_h__
#define __extfield_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for external fields
  /*! @ingroup fermstates
   *
   * Supports external field creation
   */
  class ExternalField
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~ExternalField() {}

    //! Produce the external field component
    virtual LatticeComplex operator()(int dir) const = 0;
  };

}


#endif
