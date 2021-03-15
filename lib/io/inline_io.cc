#include "chromabase.h"
#include "io/inline_io.h"
#include "meas/inline/abs_inline_measurement_factory.h"


namespace Chroma { 

  // Read an inline measurement
  void read(XMLReader& xml,
	    const std::string& path,
	    Handle< AbsInlineMeasurement >& meas_handle) 
  {
    meas_handle = readInlineMeasurement(xml, path);
  }


  // Return an inline measurement
  AbsInlineMeasurement* readInlineMeasurement(XMLReader& xml,
					      const std::string& path)
  {

    std::string measurement_name;
    try { 
      XMLReader paramtop(xml, path);
      read( paramtop, "./Name", measurement_name);
    }
    catch(const std::string& e) { 
      QDPIO::cerr << "Caught Exception Reading XML: " << e << std::endl;
      QDP_abort(1);
    }
    
    return TheInlineMeasurementFactory::Instance().createObject(measurement_name, 
								xml,
								path);
    
  }
  
}
