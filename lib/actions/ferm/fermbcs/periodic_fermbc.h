// -*- C++ -*-
// $Id: periodic_fermbc.h,v 2.2 2006-03-16 03:00:12 edwards Exp $
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
  template<class T>
  class PeriodicFermBC : public FermBC<T>
  {
  public:
    //! Only full constructor
    PeriodicFermBC() {}

    //! Copy constructor
    PeriodicFermBC(const PeriodicFermBC& a) {}

    //! Destructor is automatic
    ~PeriodicFermBC() {}

    //! Assignment
    PeriodicFermBC& operator=(const PeriodicFermBC&) {return *this;}

    //! Modify U fields in place
    /*! NOP */
    void modifyU(multi1d<LatticeColorMatrix>& u) const {}

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(T& psi) const {}
 
    //! Zero some gauge-like field in place on the masked links
    void zero(multi1d<LatticeColorMatrix>& ds_u) const {}

    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return false;}
  };

}

#endif
