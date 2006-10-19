// -*- C++ -*-
// $Id: fermact_factory_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __fermact_factory_w_h__
#define __fermact_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "wilstype_fermact_w.h"

namespace Chroma
{
  //! Wilson-like fermion factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<FermionAction<LatticeFermion,
			      multi1d<LatticeColorMatrix>,
			      multi1d<LatticeColorMatrix> >,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
						     const std::string&), StringFactoryError> >
  TheFermionActionFactory;


  //! Wilson-like fermion 4D factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct<LatticeFermion, 
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    WilsonTypeFermAct<LatticeFermion, 
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
							 const std::string&), 
		StringFactoryError> >
  TheWilsonTypeFermActFactory;


  //! Wilson-like fermion array factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<WilsonTypeFermAct5D<LatticeFermion, 
				    multi1d<LatticeColorMatrix>, 
				    multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    WilsonTypeFermAct5D<LatticeFermion, 
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
							   const std::string&), 
		StringFactoryError> >
  TheWilsonTypeFermAct5DFactory;

}


#endif
