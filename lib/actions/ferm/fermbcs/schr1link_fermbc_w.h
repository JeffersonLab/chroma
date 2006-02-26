// -*- C++ -*-
// $Id: schr1link_fermbc_w.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! @file
 * @brief 1-Link Schroedinger function boundary conditions
 */

#if 0


#ifndef __schr1link_fermbc_w_h__
#define __schr1link_fermbc_w_h__

#include "actions/ferm/fermbcs/schroedinger_fermbc_w.h"
#include "actions/gauge/gaugebcs/schr1link_gaugebc.h"

namespace Chroma
{

  //! Name and registration
  /*! \ingroup fermbc */
  namespace WilsonTypeSchr1LinkFermBCEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Concrete class for 1-link gauge action boundary conditions with Schroedinger BC
  /*! @ingroup fermbc
   *
   *  Schroedinger BC for gauge actions that have only 1 link for padding
   *  in the decay direction - e.g., Wilson-like gauge action
   */
  template<class T>
  class Schr1LinkFermBC : public SchrFermBC<T>
  {
  public:
    //! Schroedinger BC
    /*! 
     * \param gbc_    a 1-link Schroedinger GAUGE BC
     * \param theta_  twist angle for BC
     */
    Schr1LinkFermBC(const Schr1LinkGaugeBC& gbc_, const Real& theta_) : gbc(gbc_) 
      {init(theta_);}

    //! Copy constructor
    Schr1LinkFermBC(const Schr1LinkFermBC& a) : gbc(a.gbc), mask(a.mask), fld(a.fld) {}

    //! Destructor is automatic
    ~Schr1LinkFermBC() {}

    //! Assignment
    Schr1LinkFermBC& operator=(const Schr1LinkFermBC& a)
      {gbc = a.gbc; mask = a.mask; fld = a.fld; return *this;}

    //! Modify U fields in place
    void modifyU(multi1d<LatticeColorMatrix>& u) const
      {QDP_error_exit("modifyU not implemented");}

    //! Modify fermion fields in place
    void modifyF(T& psi) const
      {QDP_error_exit("modifyF not implemented");}
 
    //! Says if there are fixed links within the lattice
    bool nontrivialP() const {return true;}

  protected:
    void init(const Real& theta);
    
  private:
    // Hide default constuctor and operator=
    Schr1LinkFermBC() {}

  private:
    Schr1LinkGaugeBC  gbc;
    multi1d<LatticeBoolean> mask;
    multi1d<LatticeColorMatrix> fld;
  };

}


#endif

#endif
