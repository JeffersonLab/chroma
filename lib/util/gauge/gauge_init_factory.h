// -*- C++ -*-
// $Id: gauge_init_factory.h,v 3.1 2007-02-04 22:06:42 edwards Exp $
/*! \file
 *  \brief Factory for producing gauge initializer objects
 */

#ifndef __gauge_init_factory_h__
#define __gauge_init_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "util/gauge/gauge_init.h"

namespace Chroma
{
  //! Gauge initialization factory (foundry)
  /*! @ingroup gauge */
  typedef SingletonHolder< 
    ObjectFactory<GaugeInit,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  GaugeInit* (*)(XMLReader&,
				 const std::string&), StringFactoryError> >
  TheGaugeInitFactory;

}


#endif
