// -*- C++ -*-
/*! \file
 * \brief Inline measurement factory
 */

#ifndef __inline_measurement_factory_h__
#define __inline_measurement_factory_h__

#include "chromabase.h"
#include "singleton.h"
#include "objfactory.h"

#include "meas/inline/abs_inline_measurement.h"


namespace Chroma 
{ 

  /*! \ingroup inline */
  typedef Chroma::SingletonHolder <
    ObjectFactory<
    AbsInlineMeasurement ,
    std::string,
    TYPELIST_2(XMLReader&, const std::string&),
    AbsInlineMeasurement* (*)(XMLReader&,
			      const std::string&), 
    StringFactoryError> >
  TheInlineMeasurementFactory;

} // End namespace

#endif

