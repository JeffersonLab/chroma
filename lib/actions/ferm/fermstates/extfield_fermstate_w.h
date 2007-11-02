// -*- C++ -*-
// $Id: extfield_fermstate_w.h,v 1.4 2007-11-02 16:34:31 kostas Exp $
/*! @file
 * @brief Fermion external field state and a creator
 */

#ifndef __extfield_fermstate_h__
#define __extfield_fermstate_h__

#include "state.h"
#include "create_state.h"
#include "handle.h"
#include "actions/ferm/fermstates/extfield.h"
#include "actions/ferm/fermstates/simple_fermstate.h"


namespace Chroma
{


   /*! @ingroup fermstates */
  namespace CreateExtFieldFermStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

  //! External field state
  /*! @ingroup fermstates
   *
   * Only needs to hold an original state and a modified state
   */
  template<typename T, typename P, typename Q>
  class ExtFieldFermState : public FermState<T,P,Q>
  {
  public:
    //! Full constructor
    ExtFieldFermState(Handle< FermBC<T,P,Q> > fbc_, 
		      Handle< ExternalField >  ext_field,
		      const Q& q_) : 
      fbc(fbc_)
      {
	// Original fields
	fs = new SimpleFermState<T,P,Q>(fbc, q_);

	// Fold in U(1) fields
	Q q_u1 = q_;
	for(int mu=0; mu < q_u1.size(); ++mu)
	  q_u1[mu] *= (*ext_field)(mu);

	fs_u1 = new SimpleFermState<T,P,Q>(fbc, q_u1);
      }

    //! Destructor
    ~ExtFieldFermState() {}

    //! Return the link fields needed in constructing linear operators
    const Q& getLinks() const {return fs_u1->getLinks();}

    //! Return the original field state
    const Handle< FermState<T,P,Q> > getOriginalState() const {return fs;}

    //! Return the U(1) modified field state
    const Handle< FermState<T,P,Q> > getU1State() const {return fs_u1;}

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    ExtFieldFermState() {}  // hide default constructur
    void operator=(const ExtFieldFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    Handle< FermState<T,P,Q> > fs;
    Handle< FermState<T,P,Q> > fs_u1;
  };



  //! Create a simple ferm connection state
  /*! @ingroup fermstates
   *
   * This is a factory class for producing a connection state
   */
  template<typename T, typename P, typename Q>
  class CreateExtFieldFermState : public CreateFermState<T,P,Q>
  {
  public:
    //! Full constructor
    CreateExtFieldFermState(Handle< FermBC<T,P,Q> > fbc_,
			    Handle< ExternalField > ext_field_) : 
      fbc(fbc_), ext_field(ext_field_) {}

    //! Destructor
    ~CreateExtFieldFermState() {}
   
    //! Construct a ConnectState
    ExtFieldFermState<T,P,Q>* operator()(const Q& q) const
      {
	return new ExtFieldFermState<T,P,Q>(fbc, ext_field, q);
      }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    CreateExtFieldFermState() {}  // hide default constructur
    void operator=(const CreateExtFieldFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    Handle< ExternalField > ext_field;
  };

}


#endif
