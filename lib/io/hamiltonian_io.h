
#ifndef HAMILTONIAN_IO_H
#define HAMILTONIAN_IO_H

#include "chromabase.h"
#include "io/monomial_io.h"
#include "update/molecdyn/hamiltonian/exact_hamiltonian.h"

using namespace QDP;
using namespace Chroma; 

namespace Chroma { 

  void read(XMLReader& xml, const std::string& path, 
	    Handle< ExactLatColMatHamiltonian >& H_handle);

}


#endif 
