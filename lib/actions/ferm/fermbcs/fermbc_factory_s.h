// -*- C++ -*-
// $Id: fermbc_factory_s.h,v 1.1 2004-12-24 04:23:20 edwards Exp $
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
  typedef SingletonHolder< 
    ObjectFactory<FermBC<LatticeStaggeredFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeStaggeredFermion>* (*)(XMLReader&, const std::string&), 
		  StringFactoryError> >
  TheStaggeredTypeFermBCFactory;

}


#endif
