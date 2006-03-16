// -*- C++ -*-
// $Id: fermbc.h,v 2.2 2006-03-16 02:57:29 edwards Exp $
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
 
    //! Zero some gauge-like field in place on the masked links
    /*! 
     * This routine may be dropped in favor of zero of a template type,
     * namely the conjugate momenta.
     */
    virtual void zero(multi1d<LatticeColorMatrix>& ds_u) const = 0;

    //! Says if there are fermion non-trivial 
    virtual bool nontrivialP() const = 0;
  };

}


#endif
