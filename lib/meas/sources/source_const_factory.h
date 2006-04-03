// -*- C++ -*-
// $Id: source_const_factory.h,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Factory for producing quark prop sources
 */

#ifndef __source_const_factory_w_h__
#define __source_const_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/sources/source_construction.h"

namespace Chroma
{
  //! Propagator source factory (foundry)
  /*! @ingroup sources */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSourceConstruction<LatticePropagator>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QuarkSourceConstruction<LatticePropagator>* (*)(XMLReader&,
								  const std::string&), StringFactoryError> >
  ThePropSourceConstructionFactory;


  //! Propagator source factory (foundry)
  /*! @ingroup sources */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSourceConstruction<LatticeFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QuarkSourceConstruction<LatticeFermion>* (*)(XMLReader&,
							       const std::string&), StringFactoryError> >
  TheFermSourceConstructionFactory;

}


#endif
