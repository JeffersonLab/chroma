// -*- C++ -*-
// $Id: monomial_factory_w.h,v 1.1 2004-12-21 15:42:22 bjoo Exp $
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
    ExactUnprecWilsonTypeFermMonomial< multi1d<LatticeColorMatrix>, 
				       multi1d<LatticeColorMatrix>, 
				       LatticeFermion>, 
    std::string,
    TYPELIST_3(Handle< FermBC<LatticeFermion> >, XMLReader&, const std::string&),

      ExactUnprecWilsonTypeFermMonomial< multi1d<LatticeColorMatrix>, 
					 multi1d<LatticeColorMatrix>,
					 LatticeFermion>* (*)(Handle< FermBC<LatticeFermion> >, 
							    XMLReader&,
							    const std::string&), StringFactoryError> >
  TheExactUnprecWilsonTypeFermMonomialFactory;

  //! A factory for exact UnprecWilsonTypeFermMonomials
  typedef SingletonHolder< 
  ObjectFactory<
    ExactEvenOddPrecWilsonTypeFermMonomial< multi1d<LatticeColorMatrix>, 
				       multi1d<LatticeColorMatrix>, 
				       LatticeFermion>, 
    std::string,
    TYPELIST_3(Handle< FermBC<LatticeFermion> >, XMLReader&, const std::string&),

      ExactEvenOddPrecWilsonTypeFermMonomial< multi1d<LatticeColorMatrix>, 
					      multi1d<LatticeColorMatrix>, 
					      LatticeFermion>(*)(Handle< FermBC<LatticeFermion> >, 
							    XMLReader&,
							    const std::string&), StringFactoryError> >
  TheExactEvenOddPrecWilsonTypeFermMonomialFactory;

}; // End namespace Chroma


#endif
