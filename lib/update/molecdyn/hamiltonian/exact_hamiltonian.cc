#include "update/molecdyn/hamiltonian/exact_hamiltonian.h"
#include <typeinfo> 

namespace Chroma { 


 ExactHamiltonianParams::ExactHamiltonianParams(XMLReader& xml, const std::string& path) 
 {
   monomial_ids.resize(0);
   try { 
     XMLReader paramtop(xml, path);
     read(paramtop, "monomial_ids", monomial_ids);
   }
   catch(const std::string& e) {
     QDPIO::cout << "Caught Exception Reading XML: " << e << endl;
     QDP_abort(1);
   }
 }
  

 void read(XMLReader& xml, const std::string& path, ExactHamiltonianParams& p)
 {
   ExactHamiltonianParams tmp(xml,path);
   p = tmp;
 }

 void write(XMLWriter& xml, const std::string& path, const ExactHamiltonianParams& p)
 {
   push(xml, path);
   write(xml, "monomial_ids", p.monomial_ids);
   pop(xml);
 }
  
 void ExactHamiltonian::create(const multi1d<std::string>& monomial_ids)
 {
   // Convenience
   typedef Handle<Monomial< multi1d<LatticeColorMatrix>, 
                            multi1d<LatticeColorMatrix> > > MHandle;

   monomials.resize(0);
   if ( monomial_ids.size() > 0 ) { 

     // Resize array of handles
     monomials.resize( monomial_ids.size() );

     // Go through the list of IDs and try to bind them
     for(int i=0; i < monomial_ids.size(); i++) { 
       
       // This will hold the lookup result temporarily
       MHandle handle;

       // Lookup the ID first
       try { 
	 handle = 
	   TheNamedObjMap::Instance().getData<MHandle>(monomial_ids[i]);
       }
       catch(const std::string& e) { 
	 QDPIO::cout << "Lookup of " << monomial_ids[i] 
		     << " failed with error: " << e << endl;
	 QDP_abort(1);
       }

       // Now try to cast it to be an ExactMonomial handle
       try { 
	 QDPIO::cout << "ExactHamiltonian::create(): Trying to bind monomial with ID " << monomial_ids[i] << endl;
	    
	    
	 monomials[i] = handle.cast<ExactMon>();
       }
       catch( std::bad_cast ) {
	 QDPIO::cout << "Failed to cast monomial with ID: " << monomial_ids[i] << " to an ExactMonomial in ExactHamiltonian::create() " << endl;
	 QDP_abort(1);
       }
     }

   }
   else { 
     // Case when there are zero monomials to the Hamiltonian
     QDPIO::cout << "Attempting to construct Hamiltonian with 0 monomials." 
		 << endl;
     QDP_abort(1);
   }
 }
    



}
