// -*- C++ -*-
// $Id: schroedinger_fermbc_w.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#ifndef __schroedinger_fermbc_w_h__
#define __schroedinger_fermbc_w_h__

#include "fermbc.h"

namespace Chroma
{

  //! Abstract class for all gauge action boundary conditions with Schroedinger BC
  /*! @ingroup fermbc
   *
   *  Schroedinger BC implies periodic in dirs orthog to decay dir, and some
   *  kind of fixed BC in the decay dir.
   */
  template<class T>
  class SchrFermBC : public FermBC<T>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SchrFermBC() {}
  };

}


#endif
