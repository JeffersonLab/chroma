// -*- C++ -*-
// $Id: seqsource_factory_w.h,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Factory for producing quark prop sinks
 */

#ifndef __seqsrc_factory_h__
#define __seqsrc_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/hadron/hadron_seqsource.h"

namespace Chroma
{
  //! Sequential source factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<HadronSeqSource<LatticePropagator>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  HadronSeqSource<LatticePropagator>* (*)(XMLReader&,
							  const std::string&),
		  StringFactoryError> >
  TheWilsonHadronSeqSourceFactory;

}


#endif
