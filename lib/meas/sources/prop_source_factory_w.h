// -*- C++ -*-
// $Id: prop_source_factory_w.h,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Factory for producing quark prop sources
 */

#ifndef __prop_source_factory_w_h__
#define __prop_source_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/sources/source_construction.h"

namespace Chroma
{
  //! Propagator source factory (foundry)
  typedef SingletonHolder< 
  ObjectFactory<SourceConstruction<LatticePropagator>, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    SourceConstruction<LatticePropagator>* (*)(XMLReader&,
					       const std::string&), StringFactoryError> >
  ThePropSourceConstructionFactory;

}


#endif
