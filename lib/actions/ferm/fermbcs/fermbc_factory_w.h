// -*- C++ -*-
// $Id: fermbc_factory_w.h,v 3.0 2006-04-03 04:58:48 edwards Exp $
/*! \file
 *  \brief Fermion Boundary Condition factories
 */

#ifndef __fermbc_factory_w_h__
#define __fermbc_factory_w_h__

#include "singleton.h"
#include "objfactory.h"

#include "fermbc.h"


namespace Chroma
{
  //! FermBC factory
  /*! \ingroup fermbcs */
  typedef SingletonHolder< 
    ObjectFactory<FermBC<LatticeFermion, 
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeFermion,
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >* (*)(XMLReader&, 
							    const std::string&), 
		  StringFactoryError> >
  TheWilsonTypeFermBCFactory;

}


#endif
