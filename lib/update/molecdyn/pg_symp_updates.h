#ifndef symp_updates_h
#define symp_updates_h

#include "update/field_state.h"
#include "update/abs_hamiltonian.h"
#include "update/abs_symp_updates.h"

// A Concrete SymplecticUpdates class for pure gauge systems
// template it on GA so that I can use the relevant copy functions etc.
class PureGaugeSympUpdates : public AbsPureGaugeSympUpdates
{
public:
  // Constructor -- with an abstract Hamiltonian System reference
  // I actually want a copy in here, but of course I can't instantiate
  // Abstract classes. Hence what I will do is call the clone virtual 
  // function and drop the result into a handle
  //
  // I use the covariant return rule, and so do a dynamic cast down.
  PureGaugeSympUpdates(const AbsHamiltonian< multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >& Ham) : H(Ham.clone()) {}

  // Copy -- To be on the safe side, I'll clone the AbsHamiltonian
  //         although it should be OK just to do H(s.H) as the Handle
  //         does refcounting -- I think.
  PureGaugeSympUpdates(const PureGaugeSympUpdates& s) : H(s.H->clone()) {}

  // Clone: Fulfil contract. Allows to be copied through the base class
  PureGaugeSympUpdates* clone() const {
    return new PureGaugeSympUpdates(*this);
  }

  // The routine to get at the Hamiltonian
  const AbsHamiltonian<multi1d<LatticeColorMatrix>, 
		       multi1d<LatticeColorMatrix> >& getHam(void) const {
    return (*H);
  }

private:
  Handle< const AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			       multi1d<LatticeColorMatrix> > > H;

};

#endif
