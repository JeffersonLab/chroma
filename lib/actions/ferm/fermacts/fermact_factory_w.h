// -*- C++ -*-
// $Id: fermact_factory_w.h,v 1.1 2004-12-24 04:23:19 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __fermact_factory_w_h__
#define __fermact_factory_w_h__

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
    TYPELIST_2(XMLReader&, const std::string&),
    WilsonTypeFermAct<LatticeFermion>* (*)(XMLReader&,
					   const std::string&), StringFactoryError> >
  TheWilsonTypeFermActFactory;


  //! Wilson-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct< multi1d<LatticeFermion> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    WilsonTypeFermAct< multi1d<LatticeFermion> >* (*)(XMLReader&,
					              const std::string&), StringFactoryError> >
  TheWilsonTypeFermActArrayFactory;

  //! Wilson-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct< multi1d<LatticeFermion> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    EvenOddPrecWilsonTypeFermAct< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
											      const std::string&), StringFactoryError> >
  TheEvenOddPrecWilsonTypeFermActArrayFactory;


  //! DWF-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<EvenOddPrecDWFermActBaseArray<LatticeFermion>, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* (*)(XMLReader&,
						       const std::string&), StringFactoryError> >
  TheEvenOddPrecDWFermActBaseArrayFactory;


  //! DWF-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<UnprecDWFermActBaseArray<LatticeFermion>, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    UnprecDWFermActBaseArray<LatticeFermion>* (*)(XMLReader&,
						  const std::string&), StringFactoryError> >
  TheUnprecDWFermActBaseArrayFactory;

}


#endif
