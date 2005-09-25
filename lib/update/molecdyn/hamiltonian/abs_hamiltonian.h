// -*- C++ -*-
// $Id: abs_hamiltonian.h,v 2.0 2005-09-25 21:04:40 edwards Exp $
/*! \file
 * \brief Abstract Hamiltonian
 *
 * Abstract Hamiltonian
 */

#ifndef abs_hamiltonian_h
#define abs_hamiltonian_h

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "io/xmllog_io.h"

namespace Chroma 
{ 

  //! Abstract Hamiltonian
  /*! 
   * @ingroup hamilton
   *
   * A hamiltonian is essentially an aggregator for Monomials 
   * As such we get two kinds -- Exact Hamiltonians and inexact 
   * ones 
   */
  template<typename P, typename Q>
  class AbsHamiltonian {
  public:
    
    //! virtual descructor:
    virtual ~AbsHamiltonian() {}
    
    //! Compute dsdq for the system...
    //  The Category default goes through all the monomials
    //  and sums their contribution
    //
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const {
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      // Self description rule
      push(xml_out, "AbsHamiltonianForce");


      int num_terms = numMonomials();
      write(xml_out, "num_terms", num_terms);

      // Caller writes elem rule
      push(xml_out, "ForcesByMonomial");

      if( num_terms > 0 ) {
	push(xml_out, "elem");
	getMonomial(0).dsdq(F, s);
	pop(xml_out);

	for(int i=1; i < num_terms; i++) { 
	  push(xml_out, "elem");
	  P cur_F;
	  getMonomial(i).dsdq(cur_F, s);
	  F += cur_F;
	  pop(xml_out);
	}
      }
      pop(xml_out); // pop("ForcesByMonomial"
      pop(xml_out); // pop("AbsHamiltonianForce");
    }

    virtual void dsdq(P& F, 
		      const AbsFieldState<P,Q>& s, 
		      const multi1d<int>& monomial_list) const {
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      // Self description rule
      push(xml_out, "AbsHamiltonianForce");


      int num_terms = monomial_list.size();
      write(xml_out, "num_terms", num_terms);

      // Caller writes elem rule
      push(xml_out, "ForcesByMonomial");

      if( num_terms > 0 ) {
	push(xml_out, "elem");
	write(xml_out, "monomial", monomial_list[0]);

	getMonomial(monomial_list[0]).dsdq(F, s);
	pop(xml_out);

	for(int i=1; i < num_terms; i++) { 
	  push(xml_out, "elem");
	  write(xml_out, "monomial", monomial_list[i]);
	  P cur_F;
	  getMonomial(monomial_list[i]).dsdq(cur_F, s);
	  F += cur_F;
	  pop(xml_out);
	}
      }
      pop(xml_out); // pop("ForcesByMonomial"
      pop(xml_out); // pop("AbsHamiltonianForce");
    }
    
    //! Refresh pseudofermsions (if any)
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& s) {
      int num_terms = numMonomials();
      for(int i=0; i < num_terms; i++) { 
	getMonomial(i).refreshInternalFields(s);
      }
    }

    //! Copy Pseudofermions 
    //! This default goes through the monomials in order and 
    //  for each one copies the fields
    //
    // CAVEAT: Source and target hamiltonians have to have the same
    // number of terms. Also the types need to be the same otherwise
    // an exception will occur. If you don't like this behavior, you
    // should override it in your code.
    virtual void setInternalFields(const AbsHamiltonian<P,Q>& H_from) {
      int num_terms = numMonomials();

      // Minimal sanity checking -- Check the number of terms are the same
      int num_terms_from = H_from.numMonomials();
      
      if( num_terms_from != num_terms ) { 
	QDPIO::cerr << "Error in setInternalFields: Source Hamiltonian and Target Hamiltonian contain different number of monomials" << endl;
	QDP_abort(1);
      }

      for(int i=0; i < num_terms; i++) { 
	getMonomial(i).setInternalFields(H_from.getMonomial(i));
      }
    }

    //! Get the number of monomials.
    virtual int numMonomials(void) const =0;

  protected:
    
    //! Get hold of monomial with index i
    virtual Monomial<P,Q>& getMonomial(int i) const = 0;
    
    
    
  };
  
  
  //! Exact Abstract Hamiltonian
  /*! @ingroup hamilton
   *
   * Now define similar classes for exact algorithms.
   * These are basically the same as before but they can compute
   * energies too. Do these need to inerit?
   * Yes. Reason: We can always give it to an inexact algorithm through
   * a downcast. In that case the energy calculations would be hidden.
   */
  template<typename P, typename Q>
  class ExactAbsHamiltonian : public AbsHamiltonian<P, Q> 
  {
  public:
    
    //! virtual descructor:
    virtual ~ExactAbsHamiltonian() {}

    //! Compute dsdq for the system...
    //  The Category default goes through all the monomials
    //  and sums their contribution
    //
    //  s is the state, F is the computed force
    virtual void dsdq(P& F, const AbsFieldState<P,Q>& s) const {
      // Self Description rule
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      push(xml_out, "ExactAbsHamiltonianForce");


      int num_terms = numMonomials();
      write(xml_out, "num_terms", num_terms);
      
      // Caller writes elem rule
      push(xml_out, "ForcesByMonomial");
      if( num_terms > 0 ) {
	push(xml_out, "elem");
	getMonomial(0).dsdq(F, s);
	pop(xml_out);  // pop("elem");

	for(int i=1; i < num_terms; i++) { 
	  push(xml_out, "elem");
	  P cur_F;
	  getMonomial(i).dsdq(cur_F, s);
	  F += cur_F;
	  pop(xml_out);  // pop("elem");
	}
      }
      pop(xml_out); // Forces by Monomial
      pop(xml_out);  // ExactAbsHamiltonian
    }

    virtual void dsdq(P& F, 
		      const AbsFieldState<P,Q>& s, 
		      const multi1d<int>& monomial_list) const {
      XMLWriter& xml_out = TheXMLOutputWriter::Instance();
      // Self description rule
      push(xml_out, "ExactAbsHamiltonianForce");


      int num_terms = monomial_list.size();
      write(xml_out, "num_terms", num_terms);

      // Caller writes elem rule
      push(xml_out, "ForcesByMonomial");

      if( num_terms > 0 ) {
	push(xml_out, "elem");
	write(xml_out, "monomial", monomial_list[0]);

	getMonomial(monomial_list[0]).dsdq(F, s);
	pop(xml_out);

	for(int i=1; i < num_terms; i++) { 
	  push(xml_out, "elem");
	  write(xml_out, "monomial", monomial_list[i]);
	  P cur_F;
	  getMonomial(monomial_list[i]).dsdq(cur_F, s);
	  F += cur_F;
	  pop(xml_out);
	}
      }
      pop(xml_out); // pop("ForcesByMonomial"
      pop(xml_out); // pop("AbsHamiltonianForce");
    }


    //! Refresh pseudofermsions (if any)
    virtual void refreshInternalFields(const AbsFieldState<P,Q>& s) {
      int num_terms = numMonomials();
      for(int i=0; i < num_terms; i++) { 
	getMonomial(i).refreshInternalFields(s);
      }
    }
    
    //! Compute the energies 
    //! The total energy
    virtual void  mesE(const AbsFieldState<P,Q>& s, Double& KE, Double& PE) const 
      {

	// Self Description Rule
	XMLWriter& xml_out = TheXMLOutputWriter::Instance();
	push(xml_out, "mesE");

	KE = mesKE(s);
	PE = mesPE(s);

	pop(xml_out);

      }
    
    //! The Kinetic Energy
    virtual Double mesKE(const AbsFieldState<P,Q>& s) const 
      {
	XMLWriter& xml_out = TheXMLOutputWriter::Instance();
	push(xml_out, "mesKE");

	// Return 1/2 sum pi^2
	// may need to loop over the indices of P?
	Double KE=norm2(s.getP());

	write(xml_out, "KE", KE);
	pop(xml_out);  // pop(mesKE);
	return KE;
      }
    
    //! The Potential Energy 
    virtual Double mesPE(const AbsFieldState<P,Q>& s) const 
      {
	// Self Encapsulation Rule
	XMLWriter& xml_out = TheXMLOutputWriter::Instance();
	push(xml_out, "mesPE");
	// Cycle through all the monomials and compute their contribution
	int num_terms = numMonomials();

	write(xml_out, "num_terms", num_terms);
	multi1d<Double> PE_terms(num_terms);
	Double PE = Double(0);
      
	// Caller writes elem rule
	push(xml_out, "PEByMonomials");
	for(int i=0; i < num_terms; i++) {
	  push(xml_out, "elem");

	  Double tmp = getMonomial(i).S(s);
	  PE += tmp;
	  pop(xml_out);
	}
	pop(xml_out); // PEByMonomials
	pop(xml_out); // pop(mesPE);
	return PE;
      }
    
    //! Get the number of monomials.
    virtual int numMonomials(void) const = 0;
    
  protected:
    
    //! Get hold of monomial with index i
    virtual ExactMonomial<P,Q>& getMonomial(int i) const = 0;
    
    
  };


}; // End namespace Chroma
#endif
