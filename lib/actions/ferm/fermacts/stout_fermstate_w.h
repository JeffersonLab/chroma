// -*- C++ -*-
// $Id: stout_fermstate_w.h,v 1.2 2006-08-03 21:14:36 edwards Exp $

/*! @file 
 *  @brief Stout field state for stout links and a creator
 *
 *  Holds gauge fields at various smearing levels and caches some
 *  auxiliary things
 */

#ifndef _stout_fermstate_h
#define _stout_fermstate_h

#include "state.h"
#include "create_state.h"
#include "actions/ferm/fermacts/stout_fermstate_params.h"

namespace Chroma 
{
  /*! @ingroup fermacts */
  namespace CreateStoutFermStateEnv 
  { 
    extern const std::string name;
    extern const bool registered;
  }


  //! Stout field state
  /*! @ingroup fermacts
   *
   * Holds a stout smeared state
   */
  class StoutFermState : public FermState<LatticeFermion,
			 multi1d<LatticeColorMatrix>,
			 multi1d<LatticeColorMatrix> >
  {
  public: 
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Constructor only from a parameter structure
    StoutFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const StoutFermStateParams& p_,
		   const multi1d<LatticeColorMatrix>& u_);


    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

    //! Destructor is automagic
    ~StoutFermState() {}

    //! Return FAT Linke
    const multi1d<LatticeColorMatrix>& getLinks() const 
      {
	return smeared_links[params.n_smear];
      }

    //! Return Links at smearing level n (n=0 is thin, n=n_smear fattest)
    const multi1d<LatticeColorMatrix>& getLinks(int smear_level) const 
      {
	return smeared_links[smear_level];
      }

    //! Convenience function 
    const multi1d<LatticeColorMatrix>& getThinLinks() const {       
      return smeared_links[0];      
    }
    

    //! derivative of a force with respect to thin links. Recursive procedure
    void deriv(multi1d<LatticeColorMatrix>& F) const; 


  private:

    // Hide default constructor
    StoutFermState(){}

    // Do the smearing from level i to level i+1
    void smear_links(const multi1d<LatticeColorMatrix>& current,
		     multi1d<LatticeColorMatrix>& next);


    // Do the force recursion from level i+1, to level i
    void deriv_recurse(const multi1d<LatticeColorMatrix>&  F_plus,
		       multi1d<LatticeColorMatrix>& F_minus,
		       const int level) const;

#if 0
// Get rid of this

    // create function
    void create(Handle< FermBC<T,P,Q> > fbc_,
		const multi1d<LatticeColorMatrix>& u_,
		const multi2d<Real>& sm_fact_,
		const int n_smear_, 
		const multi1d<bool>& smear_in_this_dirP_);
#endif

    // create function
    void create(Handle< FermBC<T,P,Q> > fbc_,
		const StoutFermStateParams& p_,
		const multi1d<LatticeColorMatrix>& u_);


  private:
    Handle< FermBC<T,P,Q> >  fbc;
    multi1d< multi1d<LatticeColorMatrix> > smeared_links;
    StoutFermStateParams  params;
  }; // End class



  //! Create a stout ferm connection state
  /*! @ingroup fermacts
   *
   * This is a factory class for producing a connection state
   */
  class CreateStoutFermState : public CreateFermState<LatticeFermion,
			       multi1d<LatticeColorMatrix>,
			       multi1d<LatticeColorMatrix> >
  {
  public: 
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    CreateStoutFermState(Handle< FermBC<T,P,Q> > fbc_,
			 const StoutFermStateParams& p_) : 
      fbc(fbc_), params(p_) {}

    //! Destructor
    ~CreateStoutFermState() {}
   
    //! Construct a ConnectState
    StoutFermState* operator()(const Q& q) const
      {
	return new StoutFermState(fbc, params, q);
      }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    CreateStoutFermState() {}  // hide default constructur
    void operator=(const CreateStoutFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    StoutFermStateParams  params;
  };

} // end namespace

#endif
