// -*- C++ -*-
/*! \file
 *  \brief Fermion Boundary Condition factories
 */

#ifndef __fermbc_factory_s_h__
#define __fermbc_factory_s_h__

#include "singleton.h"
#include "objfactory.h"

#include "fermbc.h"


namespace Chroma
{
  //! FermBC factory
  /*! \ingroup fermbcs */
  typedef SingletonHolder< 
    ObjectFactory<FermBC<LatticeStaggeredFermion,
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeStaggeredFermion,
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >* (*)(XMLReader&, 
							    const std::string&), 
		  StringFactoryError> >
  TheStaggeredTypeFermBCFactory;

}


#endif
