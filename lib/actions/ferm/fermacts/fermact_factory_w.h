// -*- C++ -*-
// $Id: fermact_factory_w.h,v 1.2 2004-12-29 22:13:40 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __fermact_factory_w_h__
#define __fermact_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

namespace Chroma
{
  //! Wilson-like fermion factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<FermionAction<LatticeFermion>, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    FermionAction<LatticeFermion>* (*)(XMLReader&,
				       const std::string&), StringFactoryError> >
  TheFermionActionFactory;


  //! Wilson-like fermion 4D factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct<LatticeFermion>, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    WilsonTypeFermAct<LatticeFermion>* (*)(XMLReader&,
					   const std::string&), StringFactoryError> >
  TheWilsonTypeFermActFactory;


  //! Wilson-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct5D<LatticeFermion>, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    WilsonTypeFermAct5D<LatticeFermion>* (*)(XMLReader&,
					     const std::string&), StringFactoryError> >
  TheWilsonTypeFermAct5DFactory;

}


#endif
