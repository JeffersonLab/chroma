// -*- C++ -*-
// $Id: state.h,v 1.1 2004-01-08 03:11:47 edwards Exp $

/*! @file
 * @brief Support class for fermion actions and linear operators
 *
 * Holds things like color link fields and other info needed for linear
 * operators
 */

#ifndef __state_h__
#define __state_h__

#include "handle.h"

using namespace QDP;

//! Support class for fermion actions and linear operators
/*! @ingroup state
 *
 * Holds things like color link fields and other info needed for linear
 * operators. 
 */
class ConnectState
{
public:
  //! Return the link fields needed in constructing linear operators
  virtual const multi1d<LatticeColorMatrix>& getLinks() const = 0;

  //! Virtual destructor to help with cleanup;
  virtual ~ConnectState() {}
};


//! Proxy for support class for fermion actions and linear operators
/*! @ingroup state
 *
 * Holds things like color link fields and other info needed for linear
 * operators. 
 */
class ConnectStateProxy : public ConnectState
{
public:
  //! Initialize pointer with existing pointer
  /*! Requires that the pointer p is a return value of new */
  explicit ConnectStateProxy(const ConnectState* p=0) : state(p) {}

  //! Copy pointer (one more owner)
  ConnectStateProxy(const ConnectStateProxy& p) : state(p.state) {}

  //! Access the value to which the pointer refers
  const ConnectState& operator*() const {return state.operator*();}
  const ConnectState* operator->() const {return state.operator->();}

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const 
    {return state->getLinks();}

protected:
  //! Assignment
  /*! Could easily be supported, but not sure why to do so... */
  ConnectStateProxy& operator=(const ConnectStateProxy& p) {state = p.state; return *this;}

private:
  Handle<const ConnectState>  state;
};


//! Simple version of Connection-State 
/*! @ingroup state
 *
 * Only needs to handle a gauge field
 */
class SimpleConnectState : public ConnectState
{
public:
  //! Full constructor
  SimpleConnectState(const multi1d<LatticeColorMatrix>& u_) {u = u_;}

  //! Copy constructor
  SimpleConnectState(const SimpleConnectState& a) {u = a.u;}

  //! Destructor
  ~SimpleConnectState() {}

  //! Return the link fields needed in constructing linear operators
  const multi1d<LatticeColorMatrix>& getLinks() const {return u;}

private:
  SimpleConnectState() {}  // hide default constructur
  void operator=(const SimpleConnectState&) {} // hide =

private:
  multi1d<LatticeColorMatrix> u;
};


#endif
