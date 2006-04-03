// -*- C++ -*-
// $Id: boundcond.h,v 3.0 2006-04-03 04:58:43 edwards Exp $
/*! @file
 * @brief Base class for all boundary conditions
 */

#ifndef __boundcond_h__
#define __boundcond_h__

#include "chromabase.h"


namespace Chroma
{
  //! Base class for all boundary conditions
  /*! @ingroup actions
   */
  template<typename P, typename Q>
  class BoundCond
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~BoundCond() {}

    //! Apply the BC onto the coordinate fields in place
    virtual void modify(Q& u) const = 0;

    //! Zero momenta field (in this case the force) in place on masked sites and links
    virtual void zero(P& ds_u) const = 0;

    //! Says if there are fixed coordinates (links) within the lattice
    virtual bool nontrivialP() const = 0;
  };


}

#endif
