// -*- C++ -*-
/*! \file
 *  \brief External field factory
 */

#ifndef __extfield_factory_w_h__
#define __extfield_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "actions/ferm/fermstates/extfield.h"

namespace Chroma
{
  //! External field factory
  typedef Chroma::SingletonHolder< 
    ObjectFactory<ExternalField,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  ExternalField* (*)(XMLReader&,
				     const std::string&), StringFactoryError> >
  TheExternalFieldFactory;

}


#endif
