// -*- C++ -*-
// $Id: gaugebc_factory.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
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
  ObjectFactory<GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    GaugeBC<multi1d<LatticeColorMatrix>, 
	    multi1d<LatticeColorMatrix> >* (*)(XMLReader&, const std::string&), 
		StringFactoryError> >
  TheGaugeBCFactory;
}; // end namespace Chroma


#endif
