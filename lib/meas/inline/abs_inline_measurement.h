// -*- C++ -*-
// $Id: abs_inline_measurement.h,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 * \brief Abstract inline measurements
 */

#ifndef __abs_inline_measurement_h__
#define __abs_inline_measurement_h__

#include "chromabase.h"

namespace Chroma 
{ 

  /*! \ingroup inline */
  class AbsInlineMeasurement 
  {
  public:
    //! Virtual Destructor
    virtual ~AbsInlineMeasurement(void) {}

    //! Tell me how often I should measure this beastie
    virtual unsigned long getFrequency(void) const = 0;

    //! Do the measurement
    virtual void operator()(unsigned long update_no,
			    XMLWriter& xml_out) = 0;
  };

} // End namespace

#endif

