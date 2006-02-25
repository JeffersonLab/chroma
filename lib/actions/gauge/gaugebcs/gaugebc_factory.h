// -*- C++ -*-
// $Id: gaugebc_factory.h,v 2.1 2006-02-25 19:47:46 edwards Exp $
/*! \file
 *  \brief Gauge boundary condition factories
 */

#ifndef __gaugebc_factory_h__
#define __gaugebc_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "gaugebc.h"

namespace Chroma
{

  //! GaugeAct Factory 
  /*! @ingroup gaugebcs */
  typedef Chroma::SingletonHolder< 
  ObjectFactory<GaugeBC, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    GaugeBC* (*)(XMLReader&, const std::string&), 
		StringFactoryError> >
  TheGaugeBCFactory;
}; // end namespace Chroma


#endif
