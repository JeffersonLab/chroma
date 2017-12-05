// -*- C++ -*-
/*! \file
 * \brief Exact Hamiltonians
 */

#ifndef MG_EXACT_HAMILTONIAN_H
#define MG_EXACT_HAMILTONIAN_H

#include "chromabase.h"
#include "handle.h"

#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/hamiltonian/exact_hamiltonian.h"

#include "io/xmllog_io.h"
#include "io/monomial_io.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{



  /*! @ingroup hamilton */
  class MGExactHamiltonian : public ExactHamiltonian 
  {

  public:

    //! Construct from a list of std::string monomial_ids
    MGExactHamiltonian(const multi1d<std::string>& monomial_ids_):ExactHamiltonian( monomial_ids_) {mu=0;}
   
    //! Construct from a parameter structure
    MGExactHamiltonian(const ExactHamiltonianParams& p):ExactHamiltonian(p){mu=0}

    //! Copy constructor
    MGExactHamiltonian(const MGExactHamiltonian& H) : ExactHamiltonian(H),mu(H.mu) {}

    void set_dir(int d){mu=d;}
    int get_dir() const {return mu;}
    //! Destructor 
    ~MGExactHamiltonian(void) {}

    Double mesKE(const AbsFieldState< 
	       multi1d<LatticeColorMatrix>, 
	       multi1d<LatticeColorMatrix> > &s
	       ) const
    {
      START_CODE();

      /* Accumulate KE per site */
      LatticeDouble ke_per_site ;

      
      //ke_per_site = -Double(4);

      // Need to call the solver to compute \Xi = inv(D^2) P
      // then need to compute 1/2 Tr P \Xi  which should be real and positive

      //OLD STUFF FOLLOWS
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


  private:
    int mu ; // the direction we are evolving
    
  };

}

#endif
