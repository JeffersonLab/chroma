#ifndef abs_symp_updates_h
#define abs_symp_updates_h


#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"

// Interface definitions. Are these useful at all?
template<typename P, typename Q>
class SymplecticUpdates {
public:
  // Virtual destructor
  virtual ~SymplecticUpdates() {}

  // Virtual Copy like clone function
  virtual SymplecticUpdates* clone(void) const = 0;

  // Access but not mutate the hamiltonian
  virtual const AbsHamiltonian<P,Q>& getHam(void) const = 0;

  // Integrate the momenta a step of length eps
  virtual void leapP(AbsFieldState<P,Q>& s, Real eps) const = 0;

  // Integrate the coordinates a step of length eps
  virtual void leapQ(AbsFieldState<P,Q>& s, Real eps) const = 0;
};


// Now something a little more concrete so that we can work
class AbsPureGaugeSympUpdates
  : public SymplecticUpdates< multi1d<LatticeColorMatrix>, 
			      multi1d<LatticeColorMatrix> >
{
public:

  virtual ~AbsPureGaugeSympUpdates() {}

  virtual AbsPureGaugeSympUpdates* clone(void) const = 0;

  // Extension of the interface, so that we can write
 
  // Satisfy the interfaces

  // Get Hamiltonian 
  // Note for the symplectic updates, it does not need to be 
  // an Exact Hamiltonian, since we only need the force which they all do
  virtual const  AbsHamiltonian<multi1d<LatticeColorMatrix>, 
				multi1d<LatticeColorMatrix> >& getHam(void) const = 0;

  
  // Leap P
  virtual void leapP(AbsFieldState<multi1d<LatticeColorMatrix>, 
		                   multi1d<LatticeColorMatrix> >& s, 
		     Real eps) const
  {
   const  AbsHamiltonian<multi1d<LatticeColorMatrix>, 
                         multi1d<LatticeColorMatrix> >& H = getHam();

    multi1d<LatticeColorMatrix> dSdQ(Nd);
    
    // Compute the force into dSdQ
    H.dsdq(s, dSdQ);
    H.applyPBoundary(dSdQ);

    for(int mu = 0; mu < Nd; mu++) { 
      // p = p - dt dSdQ
      (s.getP())[mu] -= eps*dSdQ[mu];

      // Take traceless antihermitian part of projection onto SU(3)
      taproj((s.getP())[mu]);
    }
  }


  // Leap Q
  virtual void leapQ(AbsFieldState<multi1d<LatticeColorMatrix>, 
		                   multi1d<LatticeColorMatrix> >& s,
		     Real eps) const
  {
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;

    // Constant
    const multi1d<LatticeColorMatrix>& p_mom = s.getP();

    // Mutable
    multi1d<LatticeColorMatrix>& u = s.getQ();

    for(int mu = 0; mu < Nd; mu++) { 

      //  dt*p[mu]
      tmp_1 = eps*p_mom[mu];
     
      // tmp_1 = exp(dt*p[mu])  
      expmat(tmp_1, EXP_TWELFTH_ORDER);
      
      // tmp_2 = exp(dt*p[mu]) u[mu] = tmp_1 * u[mu]
      tmp_2 = tmp_1*u[mu];

      // u[mu] =  tmp_1 * u[mu] =  tmp_2 
      u[mu] = tmp_2;

      // Reunitarize u[mu]
      int numbad;
      reunit(u[mu], numbad, REUNITARIZE_ERROR);
    }

    // Do I need boundary conditions here?
    // getHam().applyQBoundary(s.getQ());
  }
};


#endif
