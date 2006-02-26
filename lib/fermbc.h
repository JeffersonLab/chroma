// -*- C++ -*-
// $Id: fermbc.h,v 2.1 2006-02-26 03:47:51 edwards Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#ifndef __fermbc_h__
#define __fermbc_h__

#include "chromabase.h"

namespace Chroma
{
  //! Base class for all gauge action boundary conditions
  /*! @ingroup fermbc
   *
   */
  template<class T>
  class FermBC
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermBC() {}

    //! Modify U fields according to the fermion BC in place
    virtual void modifyU(multi1d<LatticeColorMatrix>& u) const = 0;

    //! Modify fermion fields in place
    virtual void modifyF(T& psi) const = 0;
 
    //! Says if there are fermion non-trivial 
    virtual bool nontrivialP() const = 0;
  };

}


#endif
