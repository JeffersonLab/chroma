#include "hamiltonian_io.h"

namespace Chroma { 

 void read(XMLReader& xml, const std::string& path, 
	   Handle< ExactLatColMatHamiltonian>& H_handle)
 {

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
   ExactLatColMatHamiltonian* tmp = new ExactLatColMatHamiltonian( monomial_array );
   if( tmp == 0 ) { 
     QDPIO::cerr << "Failed to create Hamiltonian " << endl;
     QDP_abort(1);
   }

   H_handle = tmp;
 }
};
