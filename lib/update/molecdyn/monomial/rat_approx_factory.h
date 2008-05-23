// -*- C++ -*-
// $Id: rat_approx_factory.h,v 3.1 2008-05-23 21:31:34 edwards Exp $
/*! \file
 *  \brief Rational approximation factories
 */

#ifndef __rat_approx_factory_w_h__
#define __rat_approx_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "update/molecdyn/monomial/rat_approx.h"

namespace Chroma
{
  //! A factory for creating rational approximation
  /*! @ingroup monomial */
  typedef SingletonHolder< 
    ObjectFactory<RationalApprox,
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
   
		  RationalApprox* (*)(XMLReader&, const std::string&), 
		  StringFactoryError> >
  TheRationalApproxFactory;

} // End namespace Chroma


#endif
