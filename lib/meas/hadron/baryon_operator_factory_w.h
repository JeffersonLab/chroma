// -*- C++ -*-
// $Id: baryon_operator_factory_w.h,v 1.2 2006-11-22 04:17:02 juge Exp $
/*! \file
 *  \brief Factory for producing baryon operators
 */

#ifndef __baryon_operator_factory_w_h__
#define __baryon_operator_factory_w_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"
#include "meas/hadron/baryon_operator.h"

namespace Chroma
{
  //! Sequential source factory (foundry)
  /*! @ingroup hadron */
  typedef SingletonHolder< 
    ObjectFactory<BaryonOperator<LatticeFermion>, 
		  std::string,
		  TYPELIST_3(XMLReader&, const std::string&, const multi1d<LatticeColorMatrix>&),
		  BaryonOperator<LatticeFermion>* (*)(XMLReader&,
						      const std::string&,
						      const multi1d<LatticeColorMatrix>&),
		  StringFactoryError> >
  TheWilsonBaryonOperatorFactory;
}


#endif
