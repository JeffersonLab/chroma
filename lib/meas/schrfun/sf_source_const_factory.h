// -*- C++ -*-
// $Id: sf_source_const_factory.h,v 1.1 2007-08-25 04:07:41 edwards Exp $
/*! \file
 *  \brief Factory for producing quark prop sources in Schroedinger Functional
 */

#ifndef __sf_source_const_factory_w_h__
#define __sf_source_const_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/schrfun/sf_source_construction.h"

namespace Chroma
{
  //! Propagator source factory (foundry)
  /*! @ingroup schrfun */
  typedef SingletonHolder< 
    ObjectFactory<SFSourceConstruction<LatticePropagator>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  SFSourceConstruction<LatticePropagator>* (*)(XMLReader&,
							       const std::string&), StringFactoryError> >
  TheSFSourceConstructionFactory;

}


#endif
