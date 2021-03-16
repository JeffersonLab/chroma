// -*- C++ -*-
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
  typedef Chroma::SingletonHolder< 
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
  typedef Chroma::SingletonHolder< 
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
