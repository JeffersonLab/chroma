// -*- C++ -*-
// $Id: spin_insertion_factory.h,v 1.2 2006-12-04 20:38:13 edwards Exp $
/*! \file
 *  \brief Factory for producing spin insertion objects
 */

#ifndef __spin_insertion_factory_h__
#define __spin_insertion_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/hadron/spin_insertion.h"

namespace Chroma
{
  //! Spin insertion factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<SpinInsertion<LatticePropagator>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  SpinInsertion<LatticePropagator>* (*)(XMLReader&,
							const std::string&), StringFactoryError> >
  ThePropSpinInsertionFactory;


  //! Spin insertion factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<SpinInsertion<LatticeFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  SpinInsertion<LatticeFermion>* (*)(XMLReader&,
						     const std::string&), StringFactoryError> >
  TheFermSpinInsertionFactory;

}


#endif
