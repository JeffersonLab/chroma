// -*- C++ -*-
// $Id: monomial_factory_w.h,v 1.3 2004-12-23 14:55:35 bjoo Exp $
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
  //! A factory for exact UnprecWilsonTypeFermMonomials
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

}; // End namespace Chroma


#endif
