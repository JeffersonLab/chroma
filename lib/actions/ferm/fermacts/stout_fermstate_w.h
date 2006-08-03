// -*- C++ -*-
// $Id: stout_fermstate_w.h,v 1.1 2006-08-03 18:55:28 edwards Exp $

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

namespace Chroma 
{
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

    //! Explicitly specify everything. This constructor 
    //  may be made private later
    StoutFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const multi1d<LatticeColorMatrix>& u_,
		   const multi2d<Real>& sm_fact_,
		   const int n_smear_, 
		   const multi1d<bool>& smear_in_this_dirP_);

    //! Explicitly specify smearing factor tensor
    StoutFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const multi1d<LatticeColorMatrix>& u_,
		   const multi2d<Real>& sm_fact_,
		   const int n_smear_);

    //! Construct isotropic smearing in all 4 directions
    StoutFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const multi1d<LatticeColorMatrix>& u_, 
		   const Real& sm_fact_, 
		   const int   n_smear_);

    //! Construct isotopic smearing in 3 directions
    StoutFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const multi1d<LatticeColorMatrix>& u_,
		   const Real& sm_fact_, 
		   const int   n_smear_,
		   const int   j_decay);

    //! Destructor is automagic
    ~StoutFermState() {}

#if 0
#warning "RGE: no idea why this stuff is here - a copy should be a copy, not a new create"

    // Copy Constructor
    StoutFermState(const StoutFermState& s) 
      {
	START_CODE();
	create(s.getThinLinks(),
	       s.rho,
	       s.n_smear, 
	       s.smear_in_this_dirP);
	END_CODE();
      }

    // Assignment
    StoutFermState& operator=(const StoutFermState& s)
      {
	START_CODE();
	create(s.getThinLinks(),
	       s.rho,
	       s.n_smear, 
	       s.smear_in_this_dirP);

	END_CODE();
	return *this;
      }
#endif


    //! Return FAT Linke
    const multi1d<LatticeColorMatrix>& getLinks() const 
      {
	return smeared_links[n_smear];
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

    // create function
    void create(Handle< FermBC<T,P,Q> > fbc_,
		const multi1d<LatticeColorMatrix>& u_,
		const multi2d<Real>& sm_fact_,
		const int n_smear_, 
		const multi1d<bool>& smear_in_this_dirP_);

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    multi1d< multi1d<LatticeColorMatrix> > smeared_links;
    multi2d< Real > rho;
    multi1d< bool > smear_in_this_dirP; // inelegant?
    int n_smear;
  }; // End class

} // end namespace

#endif
