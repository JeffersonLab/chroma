// -*- C++ -*-
// $Id: gauge_createstate_factory.h,v 1.1 2006-09-19 18:21:38 edwards Exp $
/*! \file
 *  \brief Gauge create state factory
 */

#ifndef __gauge_createstate_factory_h__
#define __gauge_createstate_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "create_state.h"

namespace Chroma
{

  //! CreateGaugeState Factory 
  /*! @ingroup gaugestates */
  typedef SingletonHolder< 
  ObjectFactory<CreateGaugeState<multi1d<LatticeColorMatrix>, 
				 multi1d<LatticeColorMatrix> >, 
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    CreateGaugeState<multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> >* (*)(XMLReader&, 
							const std::string&), 
		StringFactoryError> >
  TheCreateGaugeStateFactory;
} // end namespace Chroma


#endif
