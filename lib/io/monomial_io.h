// -*- C++ -*-
// $Id: monomial_io.h,v 2.1 2006-03-21 19:11:28 edwards Exp $
/*! \file
 *  \brief Monomial IO
 */

#ifndef MONOMIAL_IO_H
#define MONOMIAL_IO_H

#include "chromabase.h"
#include "update/molecdyn/monomial/abs_monomial.h"


namespace Chroma 
{ 

  /*! \ingroup io */
  void read(XMLReader& xml, 
	    const std::string& path, 
	    Handle< Monomial< multi1d<LatticeColorMatrix>, 
	                      multi1d<LatticeColorMatrix> > >& mon_handle );

  //! Read a monomial from an XML reader, use a factory to create and assign the pointer to the handle...
  /*! \ingroup io */
  void read(XMLReader& xml, 
	    const std::string& path, 
	    Handle< ExactMonomial< multi1d<LatticeColorMatrix>, 
	                           multi1d<LatticeColorMatrix> > >& mon_handle );



}; // End namespace Chroma


#endif
