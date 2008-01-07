// -*- C++ -*-
// $Id: dilution_scheme_factory.h,v 1.1 2007/12/14 06:53:42 edwards Exp 
/*! \file
 *  \brief Factory for dilution schemes
 */

#ifndef __dilution_scheme_factory_h__
#define __dilution_scheme_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/hadron/dilution_scheme.h"

namespace Chroma
{
  //! Dilution operator factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<DilutionScheme<LatticeFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  DilutionScheme<LatticeFermion>* (*)(XMLReader&,
							const std::string&),
		  StringFactoryError> >
  TheFermDilutionSchemeFactory;

}


#endif
