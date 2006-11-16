// -*- C++ -*-
// $Id: monomial_io.h,v 3.1 2006-11-16 20:39:48 bjoo Exp $
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

  //! Read a named monomial from an XML reader, use factory to create and assign the pointer to a handle in the named object map of monomial handles.
  /*! \ingroup io */
  void readNamedMonomial(XMLReader& xml,
			 const std::string& path,
			 std::string& the_monomial_id);


  //! Read an array of named monomials from an XML reader. use factory to create the monomials and put them in a named object map of monomial handles.
  /*! \ingroup io */
  void readNamedMonomialArray(XMLReader& xml, 
			      const std::string& path);


}; // End namespace Chroma


#endif
