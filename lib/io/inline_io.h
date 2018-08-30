// -*- C++ -*-
/*! \file
 *  \brief Support for inline measurements
 */

#ifndef INLINE_IO_H
#define INLINE_IO_H

#include "chromabase.h"
#include "handle.h"
#include "meas/inline/abs_inline_measurement.h"


namespace Chroma { 

  // Read an inline measurement
  void read(XMLReader& xml, 
	    const std::string& path, 
	    Handle< AbsInlineMeasurement >& meas_handle);
  

  // Return an inline measurement
  AbsInlineMeasurement* readInlineMeasurement(XMLReader& xml,
					      const std::string& path);
}

#endif
