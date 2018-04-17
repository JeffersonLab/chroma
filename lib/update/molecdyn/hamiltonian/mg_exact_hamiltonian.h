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

#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include "io/xmllog_io.h"
#include "io/monomial_io.h"
#include "meas/inline/io/named_objmap.h"
#include "actions/boson/operator/adjoint_derivative.h"
#include "actions/boson/invert/invcg_adj.h"


namespace Chroma 
{



  /*! @ingroup hamilton */
  class MGExactHamiltonian : public ExactHamiltonian 
  {

  public:

    //! Construct from a list of std::string monomial_ids
    MGExactHamiltonian(const multi1d<std::string>& monomial_ids_):ExactHamiltonian( monomial_ids_) {}
   
    //! Construct from a parameter structure
    MGExactHamiltonian(const ExactHamiltonianParams& p):ExactHamiltonian(p){}

    //! Copy constructor
    MGExactHamiltonian(const MGExactHamiltonian& H) : ExactHamiltonian(H) {}

    //void set_dir(int d){mu=d;}
    //void set_rho(Real d){rho=d;}
    //int get_dir() const {return mu;}
    //! Destructor 
    ~MGExactHamiltonian(void) {}

    Double mesKE(const AbsFieldState< 
	       multi1d<LatticeColorMatrix>, 
	       multi1d<LatticeColorMatrix> > &s
	       ) const
    {
      START_CODE();
      Double KE=zero ;
      
      // First check what if we are doing MG
      if( LCMMDIntegratorSteps::theHyperPlane::Instance().active()){
	int mu = LCMMDIntegratorSteps::theHyperPlane::Instance().getDir();
	Real rho(LCMMDIntegratorSteps::theHyperPlane::Instance().getRho());
	/* Accumulate KE per site */
	//

      
	Real RsdCG = 1.0e-12 ; // hard coded for now
	int MaxCG = 10000 ; // hard coded for now
	
	// I need to define the adjoint derivative
	LatticeDouble ke_per_site ;
	// we need to pass down mu and rho
	// rho has to come in from the xml ... where is where the fun starts
	AdjointDerivative D(mu,rho,s.getQ());
	Handle< squaredAdjointDerivative> D2 = new squaredAdjointDerivative(D) ;
	//squaredAdjointDerivative D2(D) ;
	LatticeColorMatrix P = s.getP()[mu] ;
	LatticeColorMatrix X = zero ;
	SystemSolverResults_t res = InvCG_adj(*D2,P,X,RsdCG,MaxCG) ;
	//ke_per_site = -Double(4);
	
	// Need to call the solver to compute \Xi = inv(D^2) P
	// then need to compute  Tr P^\dagger \Xi  which should be real and positive
	// I use the normalization of kinnetic energy Balint has..
	// But I do not understand the -4 (see code snipet below
	/* Now add on the local Norm2 of the momenta for each link */
	//for(int mu=0; mu < Nd; mu++) { 
	//	ke_per_site[mu] = -Double(4);
	//	ke_per_site[mu] += localNorm2(s.getP()[mu]);
	//}
	//I hope this works
	KE = real(sum(trace(adj(P)*X))) ;

      }
      else{ // call the ExactHamiltonian
	
	/* Accumulate KE per site */
	multi1d<LatticeDouble> ke_per_site(Nd);
	
      
	/* Now add on the local Norm2 of the momenta for each link */
	for(int mu=0; mu < Nd; mu++) { 
	  ke_per_site[mu] = -Double(4);
	  ke_per_site[mu] += localNorm2(s.getP()[mu]);
	}
	
	/* Sum up the differences */
	for(int mu=0; mu < Nd; mu++) { 
	  KE += sum(ke_per_site[mu]);
	}
	
      }
      
      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "mesKE");
      write(xml_out, "KE", sum(KE));
      pop(xml_out);  // pop(mesKE);


      return KE;

      END_CODE();
    }


  private:
    //    int mu ; // the direction we are evolving
    //Real rho ; // the infrared cut off
  };

}

#endif
