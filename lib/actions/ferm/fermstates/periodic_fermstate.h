// -*- C++ -*-
// $Id: periodic_fermstate.h,v 1.1 2006-09-19 17:53:37 edwards Exp $

/*! @file
 * @brief Periodic ferm state and a creator
 */

#ifndef __periodic_fermstate_h__
#define __periodic_fermstate_h__

#include "state.h"
#include "create_state.h"
#include "handle.h"
#include "actions/ferm/fermbcs/periodic_fermbc.h"

namespace Chroma
{
  //! Periodic version of FermState 
  /*! @ingroup fermstates
   *
   * Only needs to hold a ferm field and ferm bc
   */
  template<typename T, typename P, typename Q>
  class PeriodicFermState : public FermState<T,P,Q>
  {
  public:
    //! Full constructor
    PeriodicFermState(const Q& q_) : 
      fbc(new PeriodicFermBC<T,P,Q>()), q(q_) {}

    //! Destructor
    ~PeriodicFermState() {}

    //! Return the link fields needed in constructing linear operators
    const Q& getLinks() const {return q;}

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    PeriodicFermState() {}  // hide default constructur
    void operator=(const PeriodicFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    Q q;
  };



  //! Create a simple ferm connection state
  /*! @ingroup fermstates
   *
   * This is a factory class for producing a connection state
   */
  template<typename T, typename P, typename Q>
  class CreatePeriodicFermState : public CreateFermState<T,P,Q>
  {
  public:
    //! Full constructor
    CreatePeriodicFermState() : fbc(new PeriodicFermBC<T,P,Q>()) {}

    //! Destructor
    ~CreatePeriodicFermState() {}
   
    //! Construct a ConnectState
    PeriodicFermState<T,P,Q>* operator()(const Q& q) const
      {
	return new PeriodicFermState<T,P,Q>(q);
      }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    void operator=(const CreatePeriodicFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
  };

}


#endif
