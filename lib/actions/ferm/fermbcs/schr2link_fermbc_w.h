// -*- C++ -*-
// $Id: schr2link_fermbc_w.h,v 2.1 2006-02-26 03:47:52 edwards Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#if 0


#ifndef __schr2link_fermbc_h__
#define __schr2link_fermbc_h__

#include "actions/ferm/fermbcs/schroedinger_fermbc_w.h"
#include "actions/gauge/gaugebcs/schr2link_gaugebc.h"

namespace Chroma
{

  //! Name and registration
  /*! \ingroup fermbc */
  namespace WilsonTypeSchr2LinkFermBCEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Concrete class for 2-link gauge action boundary conditions with Schroedinger BC
  /*! @ingroup fermbc
   *
   *  Schroedinger BC for gauge actions that have only 2 links for padding
   *  in the decay direction - e.g., Symanzik-like gauge action
   */
  template<class T>
  class Schr2LinkFermBC : public SchrFermBC<T>
  {
  public:
    //! Schroedinger BC
    /*! 
     * \param gbc_    a 2-link Schroedinger GAUGE BC
     * \param theta_  twist angle for BC
     */
    Schr2LinkFermBC(const Schr2LinkGaugeBC& gbc_, const Real& theta_) : gbc(gbc_) 
      {init(theta_);}

    //! Copy constructor
    Schr2LinkFermBC(const Schr2LinkFermBC& a) : gbc(a.gbc), mask(a.mask), fld(a.fld) {}

    //! Destructor is automatic
    ~Schr2LinkFermBC() {}

    //! Assignment
    Schr2LinkFermBC& operator=(const Schr2LinkFermBC& a)
      {gbc = a.gbc; mask = a.mask; fld = a.fld; return *this;}

//    //! Type of Schroedinger BC
//    SchrFunType getSFermBC() const {return SchrFun;}

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
    Schr2LinkFermBC() {}

  private:
    Schr2LinkGaugeBC  gbc;
    multi1d<LatticeBoolean> mask;
    multi1d<LatticeColorMatrix> fld;
  };

}


#endif


#endif
