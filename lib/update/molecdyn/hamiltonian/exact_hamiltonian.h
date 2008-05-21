// -*- C++ -*-
// $Id: exact_hamiltonian.h,v 3.5 2008-05-21 17:07:50 bjoo Exp $
/*! \file
 * \brief Exact Hamiltonians
 */

#ifndef EXACT_HAMILTONIAN_H
#define EXACT_HAMILTONIAN_H

#include "chromabase.h"
#include "handle.h"

#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"

#include "io/xmllog_io.h"
#include "io/monomial_io.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{



  //! Parameter structure for new Hamiltonian
  /*! @ingroup hamilton */
  struct ExactHamiltonianParams { 
    
    //! Constructor 
    ExactHamiltonianParams(XMLReader& xml, const std::string& path); 
    multi1d<std::string> monomial_ids; /*!< list of monomial IDs */

  };

  //! Read the parameters for the Hamiltonian
  void read(XMLReader& xml, const std::string& path, ExactHamiltonianParams& p);

  //! Write the parameters for the Hamiltonian
  void write(XMLWriter& xml, const std::string& path, const ExactHamiltonianParams& p);


  //! The Exact Hamiltonian Class - supplies internal field refreshment and energy calculations
  /*! @ingroup hamilton */
  class ExactHamiltonian : public AbsHamiltonian< multi1d<LatticeColorMatrix>, 
		      multi1d<LatticeColorMatrix> >
  {
  public:

    //! Construct from a list of string monomial_ids
    ExactHamiltonian(const multi1d<std::string>& monomial_ids_)  {
      create(monomial_ids_);
    }
   
    //! Construct from a parameter structure
    ExactHamiltonian(const ExactHamiltonianParams& p) {
      create(p.monomial_ids);
    }

    //! Copy constructor
    ExactHamiltonian(const ExactHamiltonian& H) : monomials(H.monomials) {}

    //! Destructor 
    ~ExactHamiltonian(void) {}

    //! Internal Field Refreshment 
    void refreshInternalFields(const AbsFieldState<multi1d<LatticeColorMatrix>,multi1d<LatticeColorMatrix> >& s)  
    { 
      START_CODE();
      for(int i=0; i < monomials.size(); i++) {
	monomials[i]->refreshInternalFields(s);
      }
      END_CODE();
    }
 

    Double mesKE(const AbsFieldState< 
	       multi1d<LatticeColorMatrix>, 
	       multi1d<LatticeColorMatrix> > &s
	       ) const
    {
      START_CODE();

      /* Accumulate KE per site */
      multi1d<LatticeDouble> ke_per_site(Nd);

      
      /* Now add on the local Norm2 of the momenta for each link */
      for(int mu=0; mu < Nd; mu++) { 
	ke_per_site[mu] = -Double(4);
	ke_per_site[mu] += localNorm2(s.getP()[mu]);
      }

      /* Sum up the differences */
      Double KE=zero;
      for(int mu=0; mu < Nd; mu++) { 
	KE += sum(ke_per_site[mu]);
      }

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesKE");
      write(xml_out, "KE", sum(KE));
      pop(xml_out);  // pop(mesKE);


      return KE;

      END_CODE();
    }

    //! The Potential Energy 
    Double  mesPE(const AbsFieldState< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >& s) const 
    {
      START_CODE();

      // Self Encapsulation Rule
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesPE");
      // Cycle through all the monomials and compute their contribution
      int num_terms = monomials.size();

      write(xml_out, "num_terms", num_terms);
      Double PE=zero;

      // Caller writes elem rule
      push(xml_out, "PEByMonomials");
      for(int i=0; i < num_terms; i++) 
      {
	push(xml_out, "elem");
	Double tmp;
	tmp=monomials[i]->S(s);
	PE += tmp;
	pop(xml_out); // elem
      }
      pop(xml_out); // PEByMonomials
      pop(xml_out); // pop(mesPE);
      
      END_CODE();
      return PE;
    }

  private:
    //! Convenience 
    typedef ExactMonomial< multi1d<LatticeColorMatrix>, 
			           multi1d<LatticeColorMatrix> >  ExactMon;

    //! This creates the hamiltonian. It is similar to the 
    void create(const multi1d<std::string>& monomial_ids);
    

    multi1d< Handle<ExactMon> >  monomials;

    
  };

}

#endif
