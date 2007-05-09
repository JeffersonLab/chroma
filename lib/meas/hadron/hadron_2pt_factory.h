// -*- C++ -*-
// $Id: hadron_2pt_factory.h,v 1.1 2007-05-09 17:19:44 edwards Exp $
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
    ObjectFactory<HadronCorrelator, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  HadronCorrelator* (*)(XMLReader&,
					const std::string&),
		  StringFactoryError> >
  TheHadron2PtFactory;

}


#endif
