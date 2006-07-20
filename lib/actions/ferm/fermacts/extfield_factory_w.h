// -*- C++ -*-
// $Id: extfield_factory_w.h,v 3.1 2006-07-20 20:06:52 edwards Exp $
/*! \file
 *  \brief External field factory
 */

#ifndef __extfield_factory_w_h__
#define __extfield_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "actions/ferm/fermacts/extfield.h"

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
