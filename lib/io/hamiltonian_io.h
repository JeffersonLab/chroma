// -*- C++ -*-
// $Id: hamiltonian_io.h,v 2.1 2006-03-21 19:11:28 edwards Exp $
/*! \file
 *  \brief Hamiltonion IO
 */

#ifndef HAMILTONIAN_IO_H
#define HAMILTONIAN_IO_H

#include "chromabase.h"
#include "io/monomial_io.h"
#include "update/molecdyn/hamiltonian/exact_hamiltonian.h"


namespace Chroma 
{ 

  /*! \ingroup io */
  void read(XMLReader& xml, const std::string& path, 
	    Handle< ExactLatColMatHamiltonian >& H_handle);

}


#endif 
