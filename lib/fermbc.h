// -*- C++ -*-
// $Id: fermbc.h,v 3.3 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#ifndef __fermbc_h__
#define __fermbc_h__

#include "chromabase.h"
#include "boundcond.h"

namespace Chroma
{
  //! Base class for all fermion action boundary conditions
  /*! @ingroup fermbc
   *
   */
  template<typename T, typename P, typename Q>
  class FermBC : public BoundCond<P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermBC() {}

    //! Modify fermion fields in place
    virtual void modifyF(T& psi) const = 0;

    //! Modify fermion fields in place under a subset
    virtual void modifyF(T& psi, const Subset& s) const = 0;

    //! Modify fermion fields in place
    /*! Convenience function */
    virtual void modifyF(multi1d<T>& psi) const = 0;

    //! Modify fermion fields in place under a subset
    /*! Convenience function */
    virtual void modifyF(multi1d<T>& psi, const Subset& s) const = 0;

    //! Modify U fields according to the fermion BC in place
    virtual void modify(Q& u) const = 0;

    //! Zero some gauge-like field in place on the masked links
    virtual void zero(P& ds_u) const = 0;

    //! Says if there are fermion non-trivial 
    virtual bool nontrivialP() const = 0;
  };

}


#endif
