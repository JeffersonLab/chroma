// -*- C++ -*-
// $Id: periodic_fermbc.h,v 3.1 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#ifndef __periodic_fermbc_h__
#define __periodic_fermbc_h__

#include "fermbc.h"

namespace Chroma
{

  //! Concrete class for all fermionic actions with trivial boundary conditions
  /*! @ingroup fermbcs
   *
   *  No BC
   */
  template<class T, typename P, typename Q>
  class PeriodicFermBC : public FermBC<T,P,Q>
  {
  public:
    //! Only full constructor
    PeriodicFermBC() {}

    //! Destructor is automatic
    ~PeriodicFermBC() {}

    //! Modify U fields in place
    /*! NOP */
    void modify(Q& u) const {}

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(T& psi) const {}
 
    //! Modify fermion fields in place under a subset
    /*! NOP */
    void modifyF(T& psi, const Subset& s) const {}

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(multi1d<T>& psi) const {}
    
    //! Modify fermion fields in place under a subset
    /*! NOP */
    void modifyF(multi1d<T>& psi, const Subset& s) const {}

    //! Zero some gauge-like field in place on the masked links
    void zero(P& ds_u) const {}

    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return false;}
  };

}

#endif
