// -*- C++ -*-
// $Id: fermfactory_w.h,v 1.2 2004-09-08 04:12:21 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __fermfactory_w_h__
#define __fermfactory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

namespace Chroma
{
  //! Wilson-like fermion factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct<LatticeFermion>, 
    std::string,
    TYPELIST_3(Handle< FermBC<LatticeFermion> >, XMLReader&, const std::string&),
    WilsonTypeFermAct<LatticeFermion>* (*)(Handle< FermBC<LatticeFermion> >, 
					   XMLReader&,
					   const std::string&), DefaultFactoryError> >
  TheWilsonTypeFermActFactory;


  //! Wilson-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct< multi1d<LatticeFermion> >, 
    std::string,
    TYPELIST_3(Handle< FermBC< multi1d<LatticeFermion> > >, XMLReader&, const std::string&),
    WilsonTypeFermAct< multi1d<LatticeFermion> >* (*)(Handle< FermBC< multi1d<LatticeFermion> > >, 
					              XMLReader&,
					              const std::string&), DefaultFactoryError> >
  TheWilsonTypeFermActArrayFactory;

}


#endif
