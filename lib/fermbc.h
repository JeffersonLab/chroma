// -*- C++ -*-
// $Id: fermbc.h,v 1.11 2005-01-14 20:13:04 edwards Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#ifndef __fermbc_h__
#define __fermbc_h__

#include "gaugebc.h"
#include "handle.h"

#include "actions/gauge/gaugebcs/gaugebc_simple.h"
#include "actions/gauge/gaugebcs/gaugebc_periodic.h"
#include "actions/gauge/gaugebcs/gaugebc_schroedinger.h"

namespace Chroma
{
  //! Base class for all gauge action boundary conditions
  /*! @ingroup actions
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


  //! Abstract class for all gauge action boundary conditions with Schroedinger BC
  /*! @ingroup actions
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


  //! Concrete class for 1-link gauge action boundary conditions with Schroedinger BC
  /*! @ingroup actions
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
      {init(theta);}

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


  //! Concrete class for 2-link gauge action boundary conditions with Schroedinger BC
  /*! @ingroup actions
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
      {init(theta);}

    //! Copy constructor
    Schr2LinkFermBC(const Schr2LinkFermBC& a) : gbc(a.gbc), mask(a.mask), fld(a.fld) {}

    //! Destructor is automatic
    ~Schr2LinkFermBC() {}

    //! Assignment
    Schr2LinkFermBC& operator=(const Schr2LinkFermBC& a)
      {gbc = a.gbc; mask = a.mask; fld = a.fld; return *this;}

    //! Type of Schroedinger BC
    SchrFunType getSFermBC() const {return SchrFun;}

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


  //! Concrete class for all gauge actions with trivial boundary conditions
  /*! @ingroup actions
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


  //! Concrete class for all gauge actions with simple boundary conditions
  /*! @ingroup actions
   *
   *  Simple BC, where boundary array is multiplied on the links on the
   *  edge of the lattice
   */
  template<class T>
  class SimpleFermBC : public FermBC<T>
  {
  public:
    //! Only full constructor
    /*!
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * NOTE: there is no real reason this is of type int, could be more general
     *       like Complex
     */
    SimpleFermBC(const multi1d<int>& boundary_) : 
      gbc(new PeriodicGaugeBC()), boundary(boundary_)
      {init();}
  
    //! Use both a Gauge BC and a boundary field
    /*!
     * \param gbc       A simple Gauge BC
     * \param boundary  multiply links on edge of lattice by boundary
     *
     * The type of gbc is restricted here mainly for paranoia.
     * The type could be loosened (in fact to an abstract type), but it makes
     * no sense to apply these type FermBC to a Schr1LinkGaugeBC, for example.
     *
     * Currently, no code uses this capability...
     */
    SimpleFermBC(const SimpleGaugeBC& gbc_, const multi1d<int>& boundary_) : 
      gbc(new SimpleGaugeBC(gbc_)), boundary(boundary_)
      {init();}
  
    //! Copy constructor
    SimpleFermBC(const SimpleFermBC& a) : 
      gbc(a.gbc), boundary(a.boundary), nontriv(a.nontriv) {}

    //! Destructor is automatic
    ~SimpleFermBC() {}

    //! Assignment
    SimpleFermBC& operator=(const SimpleFermBC& a)
      {gbc = a.gbc; boundary = a.boundary; nontriv = a.nontriv; return *this;}

    //! Modify U fields in place
    void modifyU(multi1d<LatticeColorMatrix>& u) const
      {
	gbc->modify(u);   // first modify the U fields

	if (nontrivialP())
	{
	  // phases for all the directions
	  for(int m = 0; m < Nd; ++m)
	  {
	    /* u[m] *=  (coord(m) == nrow(m)-1 ) ? boundary[m] : 1 */
	    u[m] *= where(Layout::latticeCoordinate(m) == Layout::lattSize()[m]-1,
			  Integer(boundary[m]), Integer(1));
	  }
	}
      }

    //! Modify fermion fields in place
    /*! NOP */
    void modifyF(T& psi) const {}
 
    //! Says if there are non-trivial BC links
    bool nontrivialP() const {return nontriv;}

  protected:
    void init()
      {
	// Check triviality of boundary
	nontriv = false;
	for(int m = 0; m < Nd; ++m)
	  if (boundary[m] != 1)
	    nontriv |= true;
      }

  private:
    // No empty constructor
    SimpleFermBC() {}

  private:
    Handle<GaugeBC> gbc;
    multi1d<int> boundary;
    bool nontriv;
  };

}


#endif
