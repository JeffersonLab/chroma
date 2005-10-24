// -*- C++ -*-
// $Id: fermbc_factory_s.h,v 2.1 2005-10-24 05:50:13 edwards Exp $
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
  /*! \ingroup fermbc */
  typedef SingletonHolder< 
    ObjectFactory<FermBC<LatticeStaggeredFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeStaggeredFermion>* (*)(XMLReader&, const std::string&), 
		  StringFactoryError> >
  TheStaggeredTypeFermBCFactory;

}


#endif
