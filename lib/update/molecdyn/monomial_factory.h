// -*- C++ -*-
// $Id: monomial_factory.h,v 1.3 2004-12-31 17:17:04 bjoo Exp $
/*! \file
 *  \brief Monomial factories
 */

#ifndef __monomial_factory_w_h__
#define __monomial_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "update/molecdyn/abs_monomial.h"



namespace Chroma
{
  //! A factory for exact non-fermionic monomials
  typedef SingletonHolder< 
  ObjectFactory<
    Monomial< multi1d<LatticeColorMatrix>, 
	      multi1d<LatticeColorMatrix> >,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
   
    Monomial< multi1d<LatticeColorMatrix>, 
	      multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
						 const std::string&), 
    StringFactoryError> >
  TheMonomialFactory;


  /*
  //! A factory for exact fermionic monomials
  typedef SingletonHolder< 
  ObjectFactory<
    ExactFermMonomial< multi1d<LatticeColorMatrix>, 
		       multi1d<LatticeColorMatrix>,
		       LatticeFermion>,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
   
    ExactFermMonomial< multi1d<LatticeColorMatrix>, 
		       multi1d<LatticeColorMatrix>,
		       LatticeFermion>* (*)(XMLReader&,
					    const std::string&), StringFactoryError> >
  TheExactFermMonomialFactory;
  */

}; // End namespace Chroma


#endif
