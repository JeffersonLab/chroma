#ifndef ferm_symp_updates_h
#define ferm_symp_updates_h

#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "update/molecdyn/abs_symp_updates.h"

// A Concrete SymplecticUpdates class for pure gauge systems
// template it on GA so that I can use the relevant copy functions etc.
class GaugeFermSympUpdates : public AbsGaugeFermSympUpdates
{
public:
  // Constructor -- with an abstract Hamiltonian System reference
  // I actually want a copy in here, but of course I can't instantiate
  // Abstract classes. Hence what I will do is call the clone virtual 
  // function and drop the result into a handle
  //
  // I use the covariant return rule, and so do a dynamic cast down.
  GaugeFermSympUpdates(const AbsFermHamiltonian< multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix>, LatticeFermion >& Ham) : H(Ham.clone()) {}

  // Copy -- To be on the safe side, I'll clone the AbsHamiltonian
  //         although it should be OK just to do H(s.H) as the Handle
  //         does refcounting -- I think.
  GaugeFermSympUpdates(const GaugeFermSympUpdates& s) : H(s.H->clone()) {}

  // Clone: Fulfil contract. Allows to be copied through the base class
  GaugeFermSympUpdates* clone() const {
    return new GaugeFermSympUpdates(*this);
  }

  // The routine to get at the Hamiltonian
  const AbsFermHamiltonian<multi1d<LatticeColorMatrix>, 
		       multi1d<LatticeColorMatrix>, LatticeFermion >& getHam(void) const {
    return (*H);
  }

private:
  Handle< const AbsFermHamiltonian<multi1d<LatticeColorMatrix>, 
                                   multi1d<LatticeColorMatrix>, 
                                   LatticeFermion > > H;
  
};

#endif
