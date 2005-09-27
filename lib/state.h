// -*- C++ -*-
// $Id: state.h,v 2.1 2005-09-27 21:16:19 bjoo Exp $

/*! @file
 * @brief Support class for fermion actions and linear operators
 *
 * Holds things like color link fields and other info needed for linear
 * operators
 */

#ifndef __state_h__
#define __state_h__

#include "handle.h"


namespace Chroma
{
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

    /*! A virtual function to get the derivative of the state.
     *  This is useful for things like fat link states, where
     *  the derivative of the state with respect to the thin 
     *  links is complicated. 
     *  The default implementation just multiplies the accumulated
     *  force by the (thin) links, which works because
     *  \dot{U} = \pi U
     * and the pi momenta get factored out 
     * this function modifies the force term 
     */
    virtual void deriv(multi1d<LatticeColorMatrix>& F) const {
      if( F.size() != Nd ) {
	throw "F has wrong size";
      }

      for(int mu=0; mu < Nd; mu++) { 
	F[mu] = getLinks()[mu]*F[mu];
      }
    }
   
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

    /*! A function to get the derivative of the state.
     *  This is useful for things like fat link states, where
     *  the derivative of the state with respect to the thin 
     *  links is complicated. 
     *  The default implementation just multiplies the accumulated
     *  force by the (thin) links, which works because
     *  \dot{U} = \pi U
     * and the pi momenta get factored out 
     * this function modifies the force term  -- This is an override
     * of a function defined in the base class
     */
    void deriv(multi1d<LatticeColorMatrix>& F) const
    {
      // Be a good little proxy and call the deriv of your encapsulated
      // state
      state->deriv(F);
    }
    
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

    // Deriv function inherited...

  private:
    SimpleConnectState() {}  // hide default constructur
    void operator=(const SimpleConnectState&) {} // hide =

  private:
    multi1d<LatticeColorMatrix> u;
  };

}


#endif
