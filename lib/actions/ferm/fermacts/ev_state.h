// -*- C++ -*-
// $Id: ev_state.h,v 1.4 2003-12-03 04:35:23 edwards Exp $
/*! @file
 * @brief Connection state holding eigenvectors
 *
 * Holds gauge fields and eigenvectors for overlap-ish thingies
 */

#ifndef __overlap_state_h__
#define __overlap_state_h__

#include "handle.h"

using namespace QDP;

//! Support class for fermion actions and linear operators
/*! @ingroup state
 *
 * Holds a gauge field and eigenvectors for some auxilliary fermion action
 *
 * The template type is that of the eigenvectors
 */
template<typename T>
class EVConnectStateBase : public ConnectState
{
public:
  //! Return the eigenvalues
  virtual const multi1d<Real>& getEigVal() const = 0;

  //! Return the eigenvectors
  virtual const multi1d<T>& getEigVec() const = 0;

  //! Return the eigenvectors
  virtual const Real& getEigValMax() const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~EVConnectStateBase() {}
};


//! Proxy for support class for fermion actions and linear operators
/*! @ingroup state
 *
 * Holds things like color link fields and other info needed for linear
 * operators. 
 */
template<typename T>
class EVConnectStateProxy : public EVConnectStateBase<T>
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  explicit EVConnectStateProxy(const EVConnectStateBase<T>* p=0) : state(p) {}

  //! Copy pointer (one more owner)
  EVConnectStateProxy(const EVConnectStateProxy& p) : state(p.state) {}

  //! Access the value to which the pointer refers
  const EVConnectStateBase<T>& operator*() const {return state.operator*();}
  const EVConnectStateBase<T>* operator->() const {return state.operator->();}

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const 
    {state->getLinks();}

  //! Return the eigenvalues
  const multi1d<Real>& getEigVal() const
    {state->getEigVal();}

  //! Return the eigenvectors
  const multi1d<T>& getEigVec() const
    {state->getEigVec();}

  //! Return the eigenvectors
  const Real& getEigValMax() const
    {state->getEigValMax();}

protected:
  //! Assignment
  /*! Could easily be supported, but not sure why to do so... */
  EVConnectStateProxy& operator=(const EVConnectStateProxy& p) {}

private:
  Handle<const EVConnectStateBase<T> >  state;
};


//! Connection-State also containing eigenvectors
/*! @ingroup fermact
 *
 * Holds a gauge field and eigenvectors for some auxilliary fermion action
 *
 * The template type is that of the eigenvectors
 */
template<typename T>
class EVConnectState : public EVConnectStateBase<T>
{
public:
  /*! 
   * Used here is a slick way of getting the underlying word type and 
   * corresponding QDP type for some lattice object
   */
//  typedef SimpleScalar< WordType<T> >  WordBase_t;
  typedef Real  WordBase_t;    // for now, use the simpler version

  //! Partial constructor
  EVConnectState(const multi1d<LatticeColorMatrix>& u_)
    {
      u = u_;   // hhmm, for now make a copy
      eigValMax = 0;
    }

  //! Full constructor
  EVConnectState(const multi1d<LatticeColorMatrix>& u_,  
		 const multi1d<WordBase_t>& val_, 
		 const multi1d<T>& vec_,
		 const WordBase_t& val_max_)
    {
      // hhmm, for now make a copy
      u = u_;
      eigVal = val_;
      eigVec = vec_;
      eigValMax = val_max_;
    }

  //! Copy constructor
  /*! 
   * Written as a deep copy, but not expected to be used 
   */
  EVConnectState(const EVConnectState& a) : u(a.u), eigVal(a.eigVal), eigVec(a.eigVec),
					    eigValMax(a.eigValMax)
    {}

  //! Destructor
  ~EVConnectState() {}

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const {return u;}

  //! Return the eigenvalues
  const multi1d<WordBase_t>& getEigVal() const {return eigVal;}

  //! Return the eigenvectors
  const multi1d<T>& getEigVec() const {return eigVec;}

  //! Return the max eigenvalues
  const WordBase_t& getEigValMax() const {return eigValMax;}

private:
  EVConnectState() {}  // hide default constructur
  void operator=(const EVConnectState&) {} // hide =

private:
  multi1d<LatticeColorMatrix> u;
  multi1d<WordBase_t>  eigVal;
  multi1d<T>           eigVec;
  WordBase_t           eigValMax;
};


#endif
