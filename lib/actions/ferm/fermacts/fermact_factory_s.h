// -*- C++ -*-
// $Id: fermact_factory_s.h,v 3.1 2006-11-17 02:17:31 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __fermact_factory_s_h__
#define __fermact_factory_s_h__

#include "singleton.h"
#include "objfactory.h"
#include "stagtype_fermact_s.h"

namespace Chroma
{
  //! Staggered-like fermion factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<FermionAction<LatticeStaggeredFermion,
			      multi1d<LatticeColorMatrix>,
			      multi1d<LatticeColorMatrix> >,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    FermionAction<LatticeStaggeredFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
						     const std::string&), StringFactoryError> >
  TheStagFermionActionFactory;


  //! Staggered-like fermion 4D factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<StaggeredTypeFermAct<LatticeStaggeredFermion, 
				     multi1d<LatticeColorMatrix>,
				     multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    StaggeredTypeFermAct<LatticeStaggeredFermion, 
			 multi1d<LatticeColorMatrix>,
			 multi1d<LatticeColorMatrix> >* (*)(XMLReader&,
							    const std::string&), 
		StringFactoryError> >
  TheStagTypeFermActFactory;

}


#endif
