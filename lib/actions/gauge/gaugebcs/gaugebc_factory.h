// -*- C++ -*-
// $Id: gaugebc_factory.h,v 2.0 2005-09-25 21:04:31 edwards Exp $
/*! \file
 *  \brief Fermion action factories
 */

#ifndef __gaugebcfactory_h__
#define __gaugebcfactory_h__

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
