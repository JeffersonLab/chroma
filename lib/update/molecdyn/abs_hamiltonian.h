#ifndef abs_hamiltonian_h
#define abs_hamiltonian_h

#include "chromabase.h"
#include "update/field_state.h"
#include "fermact.h"

/*! An abstract hamiltonian class, for inexact algorithms.
  
    Inexact in this case means energy computation is not supported,
    (in an inexact algorithm sense -- obviously it is weird to have
    a hamiltonian where you cannot compute the energy. We may need
    to think more about this)

    This serves the following purpose. It definees 
    an interface for computing the total force 
    and can refresh the momenta,


    We don't specify how the momenta is refreshed. It is "virtual".
    HMD type algorithms will porbably use gaussian noise. 
    GHMD type algorithms will mix the previous momenta some. How
       to do that will be encoded in the derived class, probably 
       through the constructor.



    For this it needs to know the types of coordinates and the momenta
    so that it can act on the right kind of state.
*/
template<typename P, typename Q>
class AbsHamiltonian {
public:
  
  //! virtual descructor:
  virtual ~AbsHamiltonian() {}

  //! clone function for virtual copy constructs
  virtual AbsHamiltonian<P,Q>* clone(void) const = 0;

  //! Compute dsdq for the system... Not specified how to actually do this
  //  s is the state, F is the computed force
  virtual void dsdq(const AbsFieldState<P,Q>& s, Q& F) const = 0;

  //! Apply any BC's to Q
  virtual void applyQBoundary(Q& q) const = 0;

  //! Apply any BC's to P
  virtual void applyPBoundary(P& p) const = 0;

};

typedef AbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > LatColMatHamiltonian;


/*! Now define similar classes for exact algorithms.
 * These are basically the same as before but they can compute
 * energies too. Do these need to inerit?
 * Yes. Reason: We can always give it to an inexact algorithm through
 * a downcast. In that case the energy calculations would be hidden.
 */

template<typename P, typename Q>
class ExactAbsHamiltonian : public AbsHamiltonian<P, Q> {
public:
  
  //! virtual descructor:
  virtual ~ExactAbsHamiltonian() {}

  //! clone function for virtual copy constructs
  virtual ExactAbsHamiltonian<P,Q>* clone(void) const = 0;

  //! Compute dsdq for the system... Not specified how to actually do this
  //  s is the state, F is the computed force
  virtual void dsdq(const AbsFieldState<P,Q>& s, Q& F) const = 0;

  //! Apply any BC's to Q
  virtual void applyQBoundary(Q& q) const = 0;

  //! Apply any BC's to P
  virtual void applyPBoundary(P& p) const = 0;


  //! Compute the energies 

  //! The total energy
  virtual void  mesE(const AbsFieldState<P,Q>& s, Double& KE, Double& PE) const {
    KE = mesKE(s);
    PE = mesPE(s);
  }

  //! The Kinetic Energy
  virtual Double mesKE(const AbsFieldState<P,Q>& s) const = 0;

  //! The Potential Energy 
  virtual Double mesPE(const AbsFieldState<P,Q>& s) const = 0;
};



template<typename P, typename Q, typename Phi>
class AbsFermHamiltonian  {
public:
  
  //! virtual descructor:
  ~AbsFermHamiltonian() {}

  //! clone function for virtual copy constructs
  virtual AbsFermHamiltonian<P,Q,Phi>* clone(void) const = 0;

  //! Get at the i-th fermion action
  virtual const FermionAction<Phi>& getFermAct(const int i) const = 0;

  //! Compute dsdq for the system... Not specified how to actually do this
  //  s is the state, F is the computed force
  virtual void dsdq(const AbsPFFieldState<P,Q,Phi>& s, Q& F) const = 0;

  virtual void applyQBoundary(Q& q) const = 0;

  //! Apply any BC's to P
  virtual void applyPBoundary(P& p) const = 0;
};


/*! Now define similar classes for exact algorithms.
 * These are basically the same as before but they can compute
 * energies too. Do these need to inerit?
 * Yes. Reason: We can always give it to an inexact algorithm through
 * a downcast. In that case the energy calculations would be hidden.
 */

template<typename P, typename Q, typename Phi>
class ExactAbsFermHamiltonian : public AbsFermHamiltonian<P, Q, Phi> {
public:
  
  //! virtual descructor:
  virtual ~ExactAbsFermHamiltonian() {}

  //! clone function for virtual copy constructs
  virtual ExactAbsFermHamiltonian<P,Q,Phi>* clone(void) const = 0;

  //! Compute dsdq for the system... Not specified how to actually do this
  //  s is the state, F is the computed force
  virtual void dsdq(const AbsPFFieldState<P,Q,Phi>& s, Q& F) const = 0;


  //! Apply any BC's to Q
  virtual void applyQBoundary(Q& q) const = 0;

  //! Apply any BC's to P
  virtual void applyPBoundary(P& p) const = 0;


  //! Get at the i-th fermion action
  virtual const FermionAction<Phi>& getFermAct(const int i) const = 0;


  //! Compute the energies 

  //! The total energy
  virtual void mesE(const AbsPFFieldState<P,Q,Phi>& s,
		    Double& KE, 
		    Double& GE, 
		    Double& FE
		   ) const {
    KE = mesKE(s);
    GE = mesGE(s);
    FE = mesFE(s);

  }

  //! The Kinetic Energy
  virtual Double mesKE(const AbsPFFieldState<P,Q,Phi>& s) const = 0;

  //! The Potential Energy 
  virtual Double mesGE(const AbsPFFieldState<P,Q,Phi>& s) const = 0;

  virtual Double mesFE(const AbsPFFieldState<P,Q,Phi>& s) const = 0;
  
};


#endif
