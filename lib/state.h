// -*- C++ -*-
// $Id: state.h,v 3.2 2006-08-26 05:40:25 edwards Exp $

/*! @file
 * @brief Support class for fermion actions and linear operators
 *
 * Holds things like color link fields and other info needed for linear
 * operators
 */

#ifndef __state_h__
#define __state_h__

#include "chromabase.h"
#include "gaugebc.h"
#include "fermbc.h"
#include "handle.h"

namespace Chroma
{
  //! Support class for fermion actions and linear operators
  /*! @ingroup state
   *
   * Holds things like color link fields and other info needed for linear
   * operators. 
   */
  template<typename P, typename Q>
  class ConnectState
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~ConnectState() {}
   
    //! Return the coordinates (link fields) needed in constructing linear operators
    virtual const Q& getLinks() const = 0;

    /*! A virtual function to get the derivative of the state.
     *  This is useful for things like fat link states, where
     *  the derivative of the state with respect to the thin 
     *  links is complicated. 
     *  The default implementation just multiplies the accumulated
     *  force by the (thin) links, which works because
     *  \f$\dot{U} = \pi U\f$
     * and the pi momenta get factored out 
     * this function modifies the force term 
     */
    virtual void deriv(P& F) const
    {
      if( F.size() != Nd ) {
	throw "F has wrong size";
      }

      P F_tmp = F;

      for(int mu=0; mu < Nd; mu++) { 
	F[mu] = getLinks()[mu]*F_tmp[mu];
      }
    }
   
    //! Return the amorphous BC object for this state
    /*! The user will supply the BC in a derived class */
    virtual const BoundCond<P,Q>& getBC() const = 0;
  };



  //! Support class for fermion actions and linear operators
  /*! @ingroup state
   *
   * Holds things like color link fields and other info needed for linear
   * operators. 
   */
  template<typename P, typename Q>
  class GaugeState : public ConnectState<P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~GaugeState() {}
   
    //! Return the gauge BC object for this state
    /*! The user will supply the BC in a derived class */
    virtual const GaugeBC<P,Q>& getBC() const = 0;
  };



  //! Support class for fermion actions and linear operators
  /*! @ingroup state
   *
   * Holds things like color link fields and other info needed for linear
   * operators. 
   */
  template<typename T, typename P, typename Q>
  class FermState : public ConnectState<P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~FermState() {}
   
    //! Return the ferm BC object for this state
    /*! The user will supply the BC in a derived class */
    virtual const FermBC<T,P,Q>& getBC() const = 0;
   
    //! Return the ferm BC object for this state
    /*! This is to help the optimized linops */
    virtual Handle< FermBC<T,P,Q> > getFermBC() const = 0;
   
  };

}


#endif
