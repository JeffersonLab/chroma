// -*- C++ -*-
// $Id: simple_fermstate.h,v 1.1 2006-09-19 17:53:37 edwards Exp $

/*! @file
 * @brief Simple ferm state and a creator
 */

#ifndef __simple_fermstate_h__
#define __simple_fermstate_h__

#include "state.h"
#include "create_state.h"
#include "handle.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"

namespace Chroma
{
  //! Simple version of FermState 
  /*! @ingroup fermstates
   *
   * Only needs to hold a ferm field and ferm bc
   */
  template<typename T, typename P, typename Q>
  class SimpleFermState : public FermState<T,P,Q>
  {
  public:
    //! Use only a SimpleFermBC with a boundary flag and a Q
    SimpleFermState(const multi1d<int>& boundary, const Q& q_) : 
      fbc(new SimpleFermBC<T,P,Q>(boundary)), q(q_)
      {
	// Apply the BC
	fbc->modify(q);
      }

    //! Full constructor
    SimpleFermState(Handle< FermBC<T,P,Q> > fbc_, const Q& q_) : 
      fbc(fbc_), q(q_) 
      {
	// Apply the BC
	fbc->modify(q);
      }

    //! Destructor
    ~SimpleFermState() {}

    //! Return the link fields needed in constructing linear operators
    const Q& getLinks() const {return q;}

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    SimpleFermState() {}  // hide default constructur
    void operator=(const SimpleFermState&) {} // hide =

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
  class CreateSimpleFermState : public CreateFermState<T,P,Q>
  {
  public:
    //! Use only a SimpleFermBC with a boundary flag
    CreateSimpleFermState(const multi1d<int>& boundary) : 
      fbc(new SimpleFermBC<T,P,Q>(boundary)) {}

    //! Full constructor
    CreateSimpleFermState(Handle< FermBC<T,P,Q> > fbc_) : fbc(fbc_) {}

    //! Destructor
    ~CreateSimpleFermState() {}
   
    //! Construct a ConnectState
    SimpleFermState<T,P,Q>* operator()(const Q& q) const
      {
	return new SimpleFermState<T,P,Q>(fbc, q);
      }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    CreateSimpleFermState() {}  // hide default constructur
    void operator=(const CreateSimpleFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
  };

}


#endif
