#include "chromabase.h"
#include "io/inline_io.h"
#include "meas/inline/abs_inline_measurement_factory.h"


namespace Chroma { 

  // Read an inline measurement
  void read(XMLReader& xml,
	    const std::string& path,
	    Handle< AbsInlineMeasurement >& meas_handle) 
  {

    std::string measurement_name;
    try { 
      XMLReader paramtop(xml, path);
      read( paramtop, "./Name", measurement_name);
    }
    catch(const std::string& e) { 
      QDPIO::cerr << "Caught Exception Reading XML: " << e << endl;
      QDP_abort(1);
    }
    
    meas_handle = TheInlineMeasurementFactory::Instance().createObject(
								      measurement_name, 
								      xml,
								      path);
    
  }
};
