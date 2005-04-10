// -*- C++ -*-
// $Id: gaugeact_factory.h,v 1.2 2005-04-10 22:32:38 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __gaugefactory_h__
#define __gaugefactory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "gaugeact.h"

namespace Chroma
{

  //! GaugeAct Factory 
  /*! @ingroup gaugeacts */
  typedef SingletonHolder< 
  ObjectFactory<GaugeAction, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    GaugeAction* (*)(XMLReader&, const std::string&), 
		StringFactoryError> >
  TheGaugeActFactory;
}; // end namespace Chroma


#endif
