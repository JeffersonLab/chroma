// -*- C++ -*-
/*! \file
 *  \brief Factory for producing hadron correlator objects
 */

#ifndef __hadron_2pt_factory_h__
#define __hadron_2pt_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/hadron/hadron_2pt.h"

namespace Chroma
{
  //! Hadron 2pt factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<HadronContract, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  HadronContract* (*)(XMLReader&,
				      const std::string&),
		  StringFactoryError> >
  TheHadronContractFactory;

}


#endif
