// -*- C++ -*-
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
  typedef Chroma::SingletonHolder< 
    ObjectFactory<LinkSmearing,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  LinkSmearing* (*)(XMLReader&,
				    const std::string&), StringFactoryError> >
  TheLinkSmearingFactory;

}


#endif
