// -*- C++ -*-
// $Id: dilution_operator_factory.h,v 1.1 2007-12-14 06:53:42 edwards Exp $
/*! \file
 *  \brief Factory for dilution objects
 */

#ifndef __dilution_operator_factory_h__
#define __dilution_operator_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/hadron/dilution_operator.h"

namespace Chroma
{
  //! Dilution operator factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<DilutionOperator<LatticeFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  DilutionOperator<LatticeFermion>* (*)(XMLReader&,
							const std::string&),
		  StringFactoryError> >
  TheFermDilutionOperatorFactory;

}


#endif
