// -*- C++ -*-
// $Id: periodic_fermbc.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#ifndef __periodic_fermbc_h__
#define __periodic_fermbc_h__

#include "fermbc.h"

namespace Chroma
{

  //! Concrete class for all fermionic actions with trivial boundary conditions
  /*! @ingroup fermbc
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
 
    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return false;}
  };

}

#endif
