// -*- C++ -*-
// $Id: monomial_factory.h,v 1.1 2004-12-29 18:38:58 bjoo Exp $
/*! \file
 *  \brief Fermion action factories
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
    ExactMonomial< multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> >,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
   
    ExactMonomial< multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
						      const std::string&), 
    StringFactoryError> >
  TheExactMonomialFactory;


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
