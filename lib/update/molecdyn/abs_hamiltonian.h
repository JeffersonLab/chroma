#ifndef abs_hamiltonian_h
#define abs_hamiltonian_h

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"


namespace Chroma { 

  /*! A hamiltonian is essentially an aggregator for Monomials 
   *  As such we get two kinds -- Exact Hamiltonians and inexact 
   *  ones */
  template<typename P, typename Q>
    class AbsHamiltonian {
    public:
    
    //! virtual descructor:
    virtual ~AbsHamiltonian() {}
    
    //! clone function for virtual copy constructs
    virtual AbsHamiltonian<P,Q>* clone(void) const = 0;
    
    
    
    //! Compute dsdq for the system...
    //  The Category default goes through all the monomials
    //  and sums their contribution
    //
    //  s is the state, F is the computed force
    virtual void dsdq(const AbsFieldState<P,Q>& s, P& F) const {
      int num_terms = numMonomials();
      
      if( num_terms > 0 ) {
	Monomial<P,Q>& first_term = getMonomial(0);
	first_term.dsdq(F, s);
	
	for(int i=1; i < num_terms; i++) { 
	  Monomial<P,Q>& current_term = getMonomial(i);
	  P cur_F;
	  current_term.dsdq(cur_F, s);
	  F += cur_F;
	}
      }
    }
    
    //! Apply any BC's to Q
    virtual void applyQBoundary(Q& q) const = 0;
    
    //! Apply any BC's to P
    virtual void applyPBoundary(P& p) const = 0;
    
    protected:
    
    //! Get hold of monomial with index i
    virtual Monomial<P,Q>& getMonomial(int i) const = 0;
    
    //! Get the number of monomials.
    virtual int numMonomials(void) const =0;
    
    
  };
  
  
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
    
    //! Apply any BC's to Q
    virtual void applyQBoundary(Q& q) const = 0;
    
    //! Apply any BC's to P
    virtual void applyPBoundary(P& p) const = 0;
    
    //! Compute dsdq for the system... Not specified how to actually do this
    //  s is the state, F is the computed force
    virtual void dsdq(const AbsFieldState<P,Q>& s, P& F) const 
    {
      int num_terms = numMonomials();
      
      if( num_terms > 0 ) {
	Monomial<P,Q>& first_term = getMonomial(0);
	first_term.dsdq(F, s);
	
	for(int i=1; i < num_terms; i++) { 
	  Monomial<P,Q>& current_term = getMonomial(i);
	  P cur_F;
	  current_term.dsdq(cur_F, s);
	  F += cur_F;
	}
      }
    }

    //! Compute the energies 
    
    //! The total energy
    virtual void  mesE(const AbsFieldState<P,Q>& s, Double& KE, Double& PE) const 
    {
      KE = mesKE(s);
      PE = mesPE(s);
    }
    
    //! The Kinetic Energy
    virtual Double mesKE(const AbsFieldState<P,Q>& s) const = 0;
    
    //! The Potential Energy 
    virtual Double mesPE(const AbsFieldState<P,Q>& s) const 
    {
	
      // Cycle through all the monomials and compute their contribution
      int num_terms = numMonomials();
      
      
      Double PE;
      if( num_terms > 0 ) {
	ExactMonomial<P,Q>& first_term = getMonomial(0);
	PE = first_term.S(s);
	
	for(int i=1; i < num_terms; i++) { 
	  ExactMonomial<P,Q>& current_term = getMonomial(i);	
	  PE += current_term.S(s);
	}
      }
    }
    
    
    protected:
    
    //! Get hold of monomial with index i
    virtual ExactMonomial<P,Q>& getMonomial(int i) const = 0;
    
    //! Get the number of monomials.
    virtual int numMonomials(void) const = 0;
    
  };

  // Potentially finger saving typedefs
  typedef AbsHamiltonian<multi1d<LatticeColorMatrix>, 
    multi1d<LatticeColorMatrix> > LatColMatHamiltonian;
  
  typedef ExactAbsHamiltonian<multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> > ExactLatColMatHamiltonian;
}; // End namespace Chroma
#endif
