// -*- C++ -*-
// $Id: gaugeact_factory.h,v 3.0 2006-04-03 04:58:54 edwards Exp $
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
  ObjectFactory<GaugeAction<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    GaugeAction<multi1d<LatticeColorMatrix>, 
		multi1d<LatticeColorMatrix> >* (*)(XMLReader&, 
						   const std::string&), 
		StringFactoryError> >
  TheGaugeActFactory;
}; // end namespace Chroma


#endif
