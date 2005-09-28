// -*- C++ -*-
// $Id: stout_state.h,v 2.2 2005-09-28 03:24:18 bjoo Exp $

/*! @file 
 *  @brief Connection State for stout links
 *
 *  Holds gauge fields at various smearing levels and caches some
 *  auxiliary things
 */

#ifndef _stout_state_h
#define _stout_state_h

#include "state.h"
#include "chromabase.h"

namespace Chroma 
{
  class StoutConnectState : public ConnectState
  {
  public: 
    
    //! Explicitly specify everything. This constructor 
    //  may be made private later
    StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
		      const multi2d<Real>& sm_fact_,
		      const int n_smear_, 
		      const multi1d<bool>& smear_in_this_dirP_);

    //! Explicitly specify smearing factor tensor
    StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
		      const multi2d<Real>& sm_fact_,
		      const int n_smear_);

    //! Construct isotropic smearing in all 4 directions
    StoutConnectState(const multi1d<LatticeColorMatrix>& u_, 
		      const Real& sm_fact_, 
		      const int   n_smear_);

    //! Construct isotopic smearing in 3 directions
    StoutConnectState(const multi1d<LatticeColorMatrix>& u_,
		      const Real& sm_fact_, 
		      const int   n_smear_,
		      const int   j_decay);

    //! Destructor is automagic
    ~StoutConnectState() {}

    // Copy Constructor
    StoutConnectState(const StoutConnectState& s) 
    {
      create(s.getThinLinks(),
	     s.rho,
	     s.n_smear, 
	     s.smear_in_this_dirP);
    }

    // Assignment
    StoutConnectState& operator=(const StoutConnectState& s)
    {
      create(s.getThinLinks(),
	     s.rho,
	     s.n_smear, 
	     s.smear_in_this_dirP);
      
      return *this;
    }


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
    StoutConnectState(){}

    // Do the smearing
    void smear_links(const multi1d<LatticeColorMatrix>& current,
		     multi1d<LatticeColorMatrix>& next);

    // create function
    void create(const multi1d<LatticeColorMatrix>& u_,
		const multi2d<Real>& sm_fact_,
		const int n_smear_, 
		const multi1d<bool>& smear_in_this_dirP_);


    multi1d< multi1d<LatticeColorMatrix> > smeared_links;
    multi2d< Real > rho;
    multi1d< bool > smear_in_this_dirP; // inelegant?
    int n_smear;
  }; // End class

} // end namespace

#endif
