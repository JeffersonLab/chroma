// -*- C++ -*-
// $Id: fermfactory_w.h,v 1.3 2004-09-18 17:12:53 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __fermfactory_w_h__
#define __fermfactory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"


namespace Chroma
{
  //! Wilson-like fermion factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct<LatticeFermion>, 
    std::string,
    TYPELIST_3(Handle< FermBC<LatticeFermion> >, XMLReader&, const std::string&),
    WilsonTypeFermAct<LatticeFermion>* (*)(Handle< FermBC<LatticeFermion> >, 
					   XMLReader&,
					   const std::string&), StringFactoryError> >
  TheWilsonTypeFermActFactory;


  //! Wilson-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct< multi1d<LatticeFermion> >, 
    std::string,
    TYPELIST_3(Handle< FermBC< multi1d<LatticeFermion> > >, XMLReader&, const std::string&),
    WilsonTypeFermAct< multi1d<LatticeFermion> >* (*)(Handle< FermBC< multi1d<LatticeFermion> > >, 
					              XMLReader&,
					              const std::string&), StringFactoryError> >
  TheWilsonTypeFermActArrayFactory;


  //! DWF-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<EvenOddPrecDWFermActBaseArray<LatticeFermion>, 
    std::string,
    TYPELIST_3(Handle< FermBC< multi1d<LatticeFermion> > >, XMLReader&, const std::string&),
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* (*)(Handle< FermBC< multi1d<LatticeFermion> > >, 
						       XMLReader&,
						       const std::string&), StringFactoryError> >
  TheEvenOddPrecDWFermActBaseArrayFactory;


  //! DWF-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<UnprecDWFermActBaseArray<LatticeFermion>, 
    std::string,
    TYPELIST_3(Handle< FermBC< multi1d<LatticeFermion> > >, XMLReader&, const std::string&),
    UnprecDWFermActBaseArray<LatticeFermion>* (*)(Handle< FermBC< multi1d<LatticeFermion> > >, 
						  XMLReader&,
						  const std::string&), StringFactoryError> >
  TheUnprecDWFermActBaseArrayFactory;

}


#endif
