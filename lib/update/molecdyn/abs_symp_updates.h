#ifndef abs_symp_updates_h
#define abs_symp_updates_h


#include "chromabase.h"
#include "update/field_state.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "util/gauge/taproj.h"
#include "util/gauge/reunit.h"
#include "util/gauge/expmat.h"

// Interface definitions. Are these useful at all?
template<typename P, typename Q, typename FS, typename H>
class AbsSymplecticUpdates {
public:
  // Virtual destructor
  virtual ~AbsSymplecticUpdates<P,Q,FS,H>() {}
  
  // Virtual copy like clone function
  virtual AbsSymplecticUpdates<P,Q,FS,H>* clone(void) const = 0;

  // access but not mutate the hamiltonian
  virtual const H& getHam(void) const = 0;

  // Integrate the momenta forward a step of length eps
  virtual void leapP(FS& s, Real eps) const = 0;

  // Integrate the coordinates forward a step of length eps
  virtual void leapQ(FS& s, Real eps) const = 0;

};

/* Specialisation to P = LatColMat, Q = LatColMat 
   at this stage we can add templated leapP and leapQ code
   which works for both fermionic and non fermionic FS and H 
*/
template<typename FS, typename H>
class AbsLatColMatSympUpdates : public 
AbsSymplecticUpdates<multi1d<LatticeColorMatrix>,
		     multi1d<LatticeColorMatrix>,
		     FS, H> {
public:

  // Virtual destructor
  virtual ~AbsLatColMatSympUpdates<FS,H>(void) {};

  // Virtual copy like clone function
  virtual AbsLatColMatSympUpdates<FS,H>* clone(void) const = 0;
  
  // Leap P
  virtual void leapP(FS& s, 
		     Real eps) const
  {
    const  H& Ham = getHam();

    multi1d<LatticeColorMatrix> dSdQ(Nd);
    
    // Compute the force into dSdQ
    Ham.dsdq(s, dSdQ);
    Ham.applyPBoundary(dSdQ);

    for(int mu = 0; mu < Nd; mu++) { 
      // p = p - dt dSdQ
      (s.getP())[mu] -= eps*dSdQ[mu];

      // Take traceless antihermitian part of projection onto SU(3)
      taproj((s.getP())[mu]);
    }
  }


  // Leap Q
  virtual void leapQ(FS& s,
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

// Now Specialise to field state and hamiltonian with no Fermions 
class AbsPureGaugeSympUpdates
  : public AbsLatColMatSympUpdates< AbsFieldState<multi1d<LatticeColorMatrix>,
					           multi1d<LatticeColorMatrix> >,
				 AbsHamiltonian<multi1d<LatticeColorMatrix>,
						multi1d<LatticeColorMatrix> > >
{
public:

  virtual ~AbsPureGaugeSympUpdates() {}

  virtual AbsPureGaugeSympUpdates* clone(void) const = 0;

  virtual const  AbsHamiltonian<multi1d<LatticeColorMatrix>, 
				multi1d<LatticeColorMatrix> >& getHam(void) const = 0;
};

class AbsGaugeFermSympUpdates
  : public AbsLatColMatSympUpdates< AbsPFFieldState<multi1d<LatticeColorMatrix>,
					           multi1d<LatticeColorMatrix>,
						    LatticeFermion >,
				 AbsFermHamiltonian<multi1d<LatticeColorMatrix>,
						    multi1d<LatticeColorMatrix>,
						    LatticeFermion> >
{
public:

  virtual ~AbsGaugeFermSympUpdates() {}

  virtual AbsGaugeFermSympUpdates* clone(void) const = 0;

  virtual const  AbsFermHamiltonian<multi1d<LatticeColorMatrix>, 
				    multi1d<LatticeColorMatrix>,
				    LatticeFermion>& getHam(void) const = 0;
};



#endif
