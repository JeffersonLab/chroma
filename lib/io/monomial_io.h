#ifndef MONOMIAL_IO_H
#define MONOMIAL_IO_H

#include "chromabase.h"
#include "update/molecdyn/abs_monomial.h"

using namespace QDP;
using namespace Chroma;

namespace Chroma { 

  void read(XMLReader& xml, 
	    const std::string& path, 
	    Handle< Monomial< multi1d<LatticeColorMatrix>, 
	                      multi1d<LatticeColorMatrix> > >& mon_handle );

  // Read a monomial from an XML reader, use a factory 
  // to create and assign the pointer to the handle...
  void read(XMLReader& xml, 
	    const std::string& path, 
	    Handle< ExactMonomial< multi1d<LatticeColorMatrix>, 
	                           multi1d<LatticeColorMatrix> > >& mon_handle );



}; // End namespace Chroma


#endif
