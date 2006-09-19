// -*- C++ -*-
// $Id: extfield_factory_w.h,v 1.1 2006-09-19 17:53:36 edwards Exp $
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
  typedef SingletonHolder< 
    ObjectFactory<ExternalField,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  ExternalField* (*)(XMLReader&,
				     const std::string&), StringFactoryError> >
  TheExternalFieldFactory;

}


#endif
