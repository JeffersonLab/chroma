#ifndef EXACT_HAMILTONIAN_H
#define EXACT_HAMILTONIAN_H

#include "chromabase.h"

#include "abs_hamiltonian.h"

using namespace QDP;
using namespace std;
using namespace Chroma;



namespace Chroma {

  class ExactLatColMatHamiltonian : public ExactAbsHamiltonian<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > {
    
  public:
    ExactLatColMatHamiltonian(XMLReader& xml, const std::string path) {
      
      XMLReader paramtop(xml,path);
      multi1d< Handle< 
	ExactMonomial<
	multi1d<LatticeColorMatrix>,
	multi1d<LatticeColorMatrix> 
	> 
	> > monomial_array;

      try { 
	read(paramtop, "./Monomials", monomial_array);
      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Error Reading Monomials " << e << endl;
	QDP_abort(1);
      }
      
      QDPIO::cout << "Read " << monomial_array.size() << " monomials" << endl;
      monomials = monomial_array;
    }

    // Constructor
    ExactLatColMatHamiltonian( multi1d< Handle< ExactMonomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > >& monomials_ ) : monomials(monomials_) {}

    // Copy Constructor
    ExactLatColMatHamiltonian( const ExactLatColMatHamiltonian& H ) : monomials(H.monomials) {}


    //! The Kinetic Energy
    Double mesKE(const AbsFieldState<multi1d<LatticeColorMatrix>,
		                     multi1d<LatticeColorMatrix> >& s) const 
    {

      // may need to loop over the indices of P?
      Double KE=Double(0);
      for(int mu=0; mu < s.getP().size(); mu++) { 
	KE += norm2((s.getP())[mu]);
      }

      QDPIO::cout << "MesKE: KE = " << KE << endl;
      return KE;
    }
    
  protected:
    ExactMonomial<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& getMonomial(int i) const {
      return *(monomials[i]);
    }

    int numMonomials(void) const {
      return monomials.size();
    }

  private:
    multi1d< Handle< ExactMonomial<multi1d<LatticeColorMatrix>,
                                   multi1d<LatticeColorMatrix> > > > monomials;
  };


};


#endif
