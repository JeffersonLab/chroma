// -*- C++ -*-
// $Id: chrono_predictor_factory.h,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 *  \brief Monomial factories
 */

#ifndef __chrono_predictor_factory_h__
#define __chrono_predictor_factory_h__

#include "singleton.h"
#include "objfactory.h"
#include "chromabase.h"

#include "update/molecdyn/predictor/chrono_predictor.h"



namespace Chroma
{
  //! A factory for exact non-fermionic monomials
  /*! @ingroup predictor */
  typedef SingletonHolder< 
  ObjectFactory< AbsChronologicalPredictor4D< LatticeFermion >,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
   
    AbsChronologicalPredictor4D<LatticeFermion>* (*)(XMLReader&,
						     const std::string&), 
    StringFactoryError> >
  The4DChronologicalPredictorFactory;

  //! A factory for exact non-fermionic monomials
  /*! @ingroup predictor */
  typedef SingletonHolder< 
  ObjectFactory< AbsChronologicalPredictor5D< LatticeFermion >,
    std::string,
    TYPELIST_3(const int, XMLReader&, const std::string&),
   
    AbsChronologicalPredictor5D<LatticeFermion>* (*)(const int,
						     XMLReader&,
						     const std::string&), 
    StringFactoryError> >
  The5DChronologicalPredictorFactory;


}; // End namespace Chroma


#endif
