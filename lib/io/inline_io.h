#ifndef INLINE_IO_H
#define INLINE_IO_H

#include "chromabase.h"
#include "handle.h"
#include "meas/inline/abs_inline_measurement.h"

using namespace QDP;
using namespace Chroma; 

namespace Chroma { 

  void read(XMLReader& xml, 
	    const std::string& path, 
	    Handle< AbsInlineMeasurement >& meas_handle);
  

};

#endif
