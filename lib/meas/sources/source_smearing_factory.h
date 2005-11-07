// -*- C++ -*-
// $Id: source_smearing_factory.h,v 2.1 2005-11-07 06:30:06 edwards Exp $
/*! \file
 *  \brief Factory for producing quark smearing objects
 */

#ifndef __source_smearing_factory_h__
#define __source_smearing_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/smear/quark_source_sink.h"

namespace Chroma
{
  //! Propagator source smearing factory (foundry)
  /*! @ingroup sources */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSourceSink<LatticePropagator>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, const multi1d<LatticeColorMatrix>&),
		  QuarkSourceSink<LatticePropagator>* (*)(XMLReader&,
							  const std::string&,
							  const multi1d<LatticeColorMatrix>&),
		  StringFactoryError> >
  ThePropSourceSmearingFactory;


  //! Propagator source smearing factory (foundry)
  /*! @ingroup sources */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSourceSink<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, const multi1d<LatticeColorMatrix>&),
		  QuarkSourceSink<LatticeFermion>* (*)(XMLReader&,
						       const std::string&,
						       const multi1d<LatticeColorMatrix>&),
		  StringFactoryError> >
  TheFermSourceSmearingFactory;

}


#endif
