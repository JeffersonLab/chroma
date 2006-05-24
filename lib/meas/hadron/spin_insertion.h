// -*- C++ -*-
// $Id: spin_insertion.h,v 1.1 2006-05-24 21:09:41 edwards Exp $
/*! @file
 * @brief Spin insertions
 */

#ifndef __spin_insertion_h__
#define __spin_insertion_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for spin insertion
  /*! @ingroup hadron
   *
   * Supports insertion into quarks
   */
  template<typename T>
  class SpinInsertion
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SpinInsertion() {}

    //! Insert into the quark
    /*!
     * \param obj      Object to insert ( Read )
     *
     * \return spin inserted object
     */
    virtual T operator()(const T& obj) const = 0;
  };

}


#endif
