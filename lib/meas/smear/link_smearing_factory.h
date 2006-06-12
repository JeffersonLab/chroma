// -*- C++ -*-
// $Id: link_smearing_factory.h,v 3.1 2006-06-12 02:13:47 edwards Exp $
/*! \file
 *  \brief Factory for producing link smearing objects
 */

#ifndef __link_smearing_factory_h__
#define __link_smearing_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/smear/link_smearing.h"

namespace Chroma
{
  //! Link smearing factory (foundry)
  /*! @ingroup smear */
  typedef SingletonHolder< 
    ObjectFactory<LinkSmearing,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  LinkSmearing* (*)(XMLReader&,
				    const std::string&), StringFactoryError> >
  TheLinkSmearingFactory;

}


#endif
