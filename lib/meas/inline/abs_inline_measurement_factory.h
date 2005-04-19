// -*- C++ -*-
// $Id: abs_inline_measurement_factory.h,v 1.2 2005-04-19 20:05:22 edwards Exp $
/*! \file
 * \brief Inline measurement factory
 */

#ifndef LCM_INLINE_MEASUREMENT_FACTORY_H
#define LCM_INLINE_MEASUREMENT_FACTORY_H

#include "chromabase.h"
#include "singleton.h"
#include "objfactory.h"

#include "meas/inline/abs_inline_measurement.h"


namespace Chroma 
{ 

  /*! \ingroup inline */
  typedef SingletonHolder <
    ObjectFactory<
    AbsInlineMeasurement ,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
   
    AbsInlineMeasurement* (*)(XMLReader&,
				       const std::string&), 
    StringFactoryError> >
  TheInlineMeasurementFactory;



}; // End namespace

#endif

