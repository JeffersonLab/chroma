// -*- C++ -*-
// $Id: gaugebc.h,v 1.3 2004-03-17 03:33:40 edwards Exp $
/*! @file
 * @brief Gauge boundary conditions
 */

#ifndef __gaugebc_h__
#define __gaugebc_h__

using namespace QDP;


//! Base class for all gauge action boundary conditions
/*! @ingroup actions
 *
 */
class GaugeBC
{
public:
  //! Virtual destructor to help with cleanup;
  virtual ~GaugeBC() {}

  //! Apply the BC onto the U fields in place
  virtual void modify(multi1d<LatticeColorMatrix>& u) const = 0;

  //! Zero the U fields in place on the masked links
  virtual void zero(multi1d<LatticeColorMatrix>& u) const = 0;

#if defined(EXPOSE_THIS_STUFF)
  // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
  //! Mask which lattice sites have fixed gauge links
  virtual const multi1d<LatticeBoolean>& lbmaskU() const = 0;

  //! Fixed gauge links on only the lbmaskU() sites
  virtual const multi1d<LatticeColorMatrix>& lFldU() const = 0;
#endif

  //! Says if there are fixed links within the lattice
  virtual bool nontrivialP() const = 0;
};


//! Schroedinger Functional type Boundary Conditions
enum SchrFunType_t {
  SF_NONE = 0, 
  SF_TRIVIAL = 1,
  SF_NONPERT = 2,
  SF_COUPLING = 3,
  SF_CHROMOMAG = 4,
  SF_DIRICHLET = 10,
};



//! Abstract class for all gauge action boundary conditions with Schroedinger BC
/*! @ingroup actions
 *
 *  Schroedinger BC implies periodic in dirs orthog to decay dir, and some
 *  kind of fixed BC in the decay dir.
 */
class SchrGaugeBC : public GaugeBC
{
public:
  //! Virtual destructor to help with cleanup;
  virtual ~SchrGaugeBC() {}

#if defined(EXPOSE_THIS_STUFF)
  //! Type of Schroedinger BC
  virtual SchrFunType_t getSFBC() const = 0;
#endif

  //! Decay direction
  virtual int getDir() const = 0;
};


//! Concrete class for 1-link gauge action boundary conditions with Schroedinger BC
/*! @ingroup actions
 *
 *  Schroedinger BC for gauge actions that have only 1 link for padding
 *  in the decay direction - e.g., Wilson-like gauge action
 */
class Schr1LinkGaugeBC : public SchrGaugeBC
{
public:
  //! Schroedinger BC
  /*! 
   * \param SchrFun_      type of Schroedinger BC
   * \param SchrPhiMult_  factor to rescale fixed field
   */
  Schr1LinkGaugeBC(SchrFunType_t SchrFun_, const Real& SchrPhiMult_);

  //! Copy constructor
  Schr1LinkGaugeBC(const Schr1LinkGaugeBC& a) : SchrFun(a.SchrFun), 
    decay_dir(a.decay_dir), mask(a.mask), fld(a.fld) {}

  //! Destructor is automatic
  ~Schr1LinkGaugeBC() {}

  //! Assignment
  Schr1LinkGaugeBC& operator=(const Schr1LinkGaugeBC& a)
    {
      SchrFun = a.SchrFun; decay_dir = a.decay_dir; 
      mask = a.mask; fld = a.fld;
      return *this;
    }

#if defined(EXPOSE_THIS_STUFF)
  //! Type of Schroedinger BC
  SchrFunType_t getSFBC() const {return SchrFun;}
#endif

  //! Modify U fields in place
  void modify(multi1d<LatticeColorMatrix>& u) const
    {QDP_error_exit("modify not implemented");}

  //! Zero the U fields in place on the masked links
  void zero(multi1d<LatticeColorMatrix>& u) const
    {QDP_error_exit("zero not implemented");}

#if defined(EXPOSE_THIS_STUFF)
  // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
  //! Mask which lattice sites have fixed gauge links
  const multi1d<LatticeBoolean>& lbmaskU() const {return mask;}

  //! Fixed gauge links on only the lbmaskU() sites
  const multi1d<LatticeColorMatrix>& lFldU() const {return fld;}
#endif

  //! Says if there are fixed links within the lattice
  bool nontrivialP() const {return true;}

  //! Decay direction
  int getDir() const {return decay_dir;}

private:
  // Hide default constuctor
  Schr1LinkGaugeBC() {}

private:
  SchrFunType_t SchrFun;
  int decay_dir;
  multi1d<LatticeBoolean> mask;
  multi1d<LatticeColorMatrix> fld;
};


//! Concrete class for 2-link gauge action boundary conditions with Schroedinger BC
/*! @ingroup actions
 *
 *  Schroedinger BC for gauge actions that have only 2 links for padding
 *  in the decay direction - e.g., Symanzik-like gauge action
 */
class Schr2LinkGaugeBC : public SchrGaugeBC
{
public:
  //! Schroedinger BC
  /*! 
   * \param SchrFun_      type of Schroedinger BC
   * \param SchrPhiMult_  factor to rescale fixed field
   */
  Schr2LinkGaugeBC(SchrFunType_t SchrFun_, const Real& SchrPhiMult_);

  //! Copy constructor
  Schr2LinkGaugeBC(const Schr2LinkGaugeBC& a) : SchrFun(a.SchrFun), 
    decay_dir(a.decay_dir), mask(a.mask), fld(a.fld) {}

  //! Destructor is automatic
  ~Schr2LinkGaugeBC() {}

  //! Assignment
  Schr2LinkGaugeBC& operator=(const Schr2LinkGaugeBC& a)
    {
      SchrFun = a.SchrFun; decay_dir = a.decay_dir; 
      mask = a.mask; fld = a.fld;
      return *this;
    }

#if defined(EXPOSE_THIS_STUFF)
  //! Type of Schroedinger BC
  SchrFunType_t getSFBC() const {return SchrFun;}
#endif

  //! Modify U fields in place
  void modify(multi1d<LatticeColorMatrix>& u) const
    {QDP_error_exit("modify not implemented");}

  //! Zero the U fields in place on the masked links
  void zero(multi1d<LatticeColorMatrix>& u) const
    {QDP_error_exit("zero not implemented");}

#if defined(EXPOSE_THIS_STUFF)
  // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
  //! Mask which lattice sites have fixed gauge links
  const multi1d<LatticeBoolean>& lbmaskU() const {return mask;}

  //! Fixed gauge links on only the lbmaskU() sites
  const multi1d<LatticeColorMatrix>& lFldU() const {return fld;}
#endif

  //! Says if there are fixed links within the lattice
  bool nontrivialP() const {return true;}

  //! Decay direction
  int getDir() const {return decay_dir;}

private:
  // Hide default constuctor
  Schr2LinkGaugeBC() {}

private:
  SchrFunType_t SchrFun;
  int decay_dir;
  multi1d<LatticeBoolean> mask;
  multi1d<LatticeColorMatrix> fld;
};


//! Concrete class for gauge actions with periodic boundary conditions
/*! @ingroup actions
 *
 *  No gauge BC - e.g., only periodic
 */
class PeriodicGaugeBC : public GaugeBC
{
public:
  //! Only empty constructor
  PeriodicGaugeBC() {}

  //! Copy constructor
  PeriodicGaugeBC(const PeriodicGaugeBC& a) {}

  //! Destructor is automatic
  ~PeriodicGaugeBC() {}

  //! Assignment
  PeriodicGaugeBC& operator=(const PeriodicGaugeBC&) {return *this;}

  //! Modify U fields in place
  /*! NOP */
  void modify(multi1d<LatticeColorMatrix>& u) const {}

  //! Zero the U fields in place on the masked links
  void zero(multi1d<LatticeColorMatrix>& u) const {}

#if defined(EXPOSE_THIS_STUFF)
  // NOT SURE THIS STUFF IS ABSOLUTELY REQUIRED - TRY TO AVOID EXPOSING THIS
  //! Mask which lattice sites have fixed gauge links
  const multi1d<LatticeBoolean>& lbmaskU() const {return mask;}

  //! Fixed gauge links on only the lbmaskU() sites
  const multi1d<LatticeColorMatrix>& lFldU() const {return fld;}
#endif

  //! Says if there are non-trivial BC links
  bool nontrivialP() const {return false;}

private:
#if defined(EXPOSE_THIS_STUFF)
  multi1d<LatticeBoolean> mask;     // is empty
  multi1d<LatticeColorMatrix> fld;  // is empty
#endif
};


//! Concrete class for gauge actions with simple boundary conditions
/*! @ingroup actions
 *
 *  Simple BC, where boundary array is multiplied on the links on the
 *  edge of the lattice
 */
class SimpleGaugeBC : public GaugeBC
{
public:
  //! Only full constructor
  /*!
   * \param boundary  multiply links on edge of lattice by boundary
   *
   * NOTE: the type of boundary can be generalized
   */
  SimpleGaugeBC(const multi1d<Complex>& boundary_) : boundary(boundary_) {}

  //! Copy constructor
  SimpleGaugeBC(const SimpleGaugeBC& a) : boundary(a.boundary) {}

  //! Destructor is automatic
  ~SimpleGaugeBC() {}

  //! Assignment
  SimpleGaugeBC& operator=(const SimpleGaugeBC& a)
    {boundary = a.boundary; return *this;}

  //! Modify U fields in place
  void modify(multi1d<LatticeColorMatrix>& u) const
    {
      // phases for all the directions
      for(int m = 0; m < Nd; ++m)
      {
	/* u[m] *=  (coord(m) == nrow(m)-1 ) ? boundary[m] : 1 */
	u[m] *= where(Layout::latticeCoordinate(m) == Layout::lattSize()[m]-1,
		      boundary[m], Complex(1));
      }
    }

  //! Zero the U fields in place on the masked links
  void zero(multi1d<LatticeColorMatrix>& u) const {}

  //! Says if there are non-trivial BC links
  /*! 
   * This is a simple implementation - always do the work. 
   * Could be improved by checking boundary
   */
  bool nontrivialP() const {return true;}

private:
  // Hide empty constructor
  SimpleGaugeBC() {}

private:
  multi1d<Complex> boundary;
};


#endif
