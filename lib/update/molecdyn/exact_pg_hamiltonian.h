#ifndef hamiltonian_h
#define hamiltonian_h

#include "chromabase.h"
#include "gaugeact.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"


//! An exact pure gauge hamiltonian. 
//! We template on the GaugeAction. If the template parameter does not
//  conform to at least the base gauge action the compiler will scream
template<typename GA>
class ExactPureGaugeHamiltonian 
  : public ExactAbsHamiltonian<multi1d<LatticeColorMatrix>,
			       multi1d<LatticeColorMatrix> >
{
public:
  //! Constructor -- takes the action as a parameter
  ExactPureGaugeHamiltonian(const GA& S_g_) : S_g(S_g_) {}

  //! Copy Constructor
  ExactPureGaugeHamiltonian(const ExactPureGaugeHamiltonian<GA>& H) : S_g(H.S_g) {}


  // Destructor
  ~ExactPureGaugeHamiltonian() {}

  // Satisfy interface:
  //! Clone Function for virtual copy
  virtual ExactPureGaugeHamiltonian* clone() const { 
    return new ExactPureGaugeHamiltonian(*this);
  }

  // Apply boundaries to Q
  virtual void applyQBoundary(multi1d<LatticeColorMatrix>& q) const
  {
    S_g.getGaugeBC().modify(q);
  }

  // Apply boundaries to P
  virtual void applyPBoundary(multi1d<LatticeColorMatrix>& p) const
  {
    S_g.getGaugeBC().zero(p);
  }

  virtual void dsdq(const AbsFieldState<multi1d<LatticeColorMatrix>,
		                        multi1d<LatticeColorMatrix> >& s,  
		    multi1d<LatticeColorMatrix>& F) const {

    // Get U from the state
    // Call CreateState to return a copy with BC's applied

    // Call dsdu from the gauge action
    S_g.dsdu(F, s.getQ());

    // ConnectState now gets freed
  }

  //! Measure the kinetic energy
  //! I am not going to normalise here, as it lets me
  //! define acceptReject in the HMC without having to know
  //! the stupid normalisations
  virtual Double mesKE(AbsFieldState<multi1d<LatticeColorMatrix>, 
		       multi1d<LatticeColorMatrix> >& s) const
  {
    // Extract momenta from the state
    // Once momenta have been refreshed, they are zeroed
    // on the gauge boundaries as requires so no need to mess
    // with gauge boundaries...
    multi1d<LatticeColorMatrix>& mom = s.getP();

    // Square it up
    Double p_mom_sq = Double(0);


    // OK Here I am not sure about Sets/Subsets. GaugeAction has a 
    // GetSet() function but it is not a "subset" so do I just play 
    // on the whole lattice? I may need to insert some subsetting code
    // in here.
    // But for now just do the do...

    // This bit here is stolen from SZIN's MesMom()
    for(int mu = 0; mu < Nd; mu++) { 
      p_mom_sq += norm2(mom[mu]);
    }
    
    return p_mom_sq;
  }

  //! Measure the potential energy (per link)
  virtual Double mesPE(AbsFieldState<multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >& s) const
  {

    // Forward the action calculation from S_g
    // There are too many S's in this line.

    // The S() function of S_g applies boundaries
    Double PE = S_g.S(s.getQ());
    return PE;
  }


private:
  // Make the constructor inaccessible
  ExactPureGaugeHamiltonian() {};

  GA S_g;          // Copy from input
};
#endif
