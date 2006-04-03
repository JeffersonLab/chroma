// -*- C++ -*-
// $Id: create_state.h,v 3.0 2006-04-03 04:58:43 edwards Exp $

/*! @file
 * @brief Create a connection state
 *
 * This is a factory class for producing a connection state
 */

#ifndef __create_state_h__
#define __create_state_h__

#include "state.h"
#include "gaugebc.h"
#include "fermbc.h"
#include "handle.h"

namespace Chroma
{
  //! Create a connection state
  /*! @ingroup state
   *
   * This is a factory class for producing a connection state
   */
  template<typename P, typename Q>
  class CreateState
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~CreateState() {}
   
    //! Construct a ConnectState
    virtual ConnectState<P,Q>* operator()(const Q& q) const = 0;

    //! Return the amorphous BC object for this state
    /*! The user will supply the BC in a derived class */
    virtual const BoundCond<P,Q>& getBC() const = 0;
  };


  //! Create a gauge connection state
  /*! @ingroup state
   *
   * This is a factory class for producing a connection state
   */
  template<typename P, typename Q>
  class CreateGaugeState : public CreateState<P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~CreateGaugeState() {}
   
    //! Construct a ConnectState
    virtual GaugeState<P,Q>* operator()(const Q& q) const = 0;

    //! Return the gauge BC object for this state
    /*! The user will supply the BC in a derived class */
    virtual const GaugeBC<P,Q>& getBC() const = 0;
  };



  //! Create a fermion connection state
  /*! @ingroup state
   *
   * This is a factory class for producing a connection state
   */
  template<typename T, typename P, typename Q>
  class CreateFermState : public CreateState<P,Q>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~CreateFermState() {}
   
    //! Construct a ConnectState
    virtual FermState<T,P,Q>* operator()(const Q& q) const = 0;

    //! Return the gauge BC object for this state
    /*! The user will supply the BC in a derived class */
    virtual const FermBC<T,P,Q>& getBC() const = 0;
   
    //! Return the gauge BC object for this state
    /*! This is to help the optimized linops */
    virtual Handle< FermBC<T,P,Q> > getFermBC() const = 0;
  };


}


#endif
