// -*- C++ -*-
// $Id: quark_smearing_factory.h,v 3.2 2006-11-17 02:17:32 edwards Exp $
/*! \file
 *  \brief Factory for producing quark smearing objects
 */

#ifndef __quark_smearing_factory_h__
#define __quark_smearing_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/smear/quark_smearing.h"

namespace Chroma
{
  //! Quark smearing factory (foundry)
  /*! \ingroup smear */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSmearing<LatticePropagator>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QuarkSmearing<LatticePropagator>* (*)(XMLReader&,
							const std::string&), StringFactoryError> >
  ThePropSmearingFactory;


  //! Quark smearing factory (foundry)
  /*! \ingroup smear */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSmearing<LatticeStaggeredPropagator>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QuarkSmearing<LatticeStaggeredPropagator>* (*)(XMLReader&,
								 const std::string&), StringFactoryError> >
  TheStagPropSmearingFactory;


  //! Quark smearing factory (foundry)
  /*! \ingroup smear */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSmearing<LatticeFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QuarkSmearing<LatticeFermion>* (*)(XMLReader&,
						     const std::string&), StringFactoryError> >
  TheFermSmearingFactory;


  //! Quark smearing factory (foundry)
  /*! \ingroup smear */
  typedef SingletonHolder< 
    ObjectFactory<QuarkSmearing<LatticeColorVector>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QuarkSmearing<LatticeColorVector>* (*)(XMLReader&,
							 const std::string&), StringFactoryError> >
  TheColorVecSmearingFactory;

}


#endif
