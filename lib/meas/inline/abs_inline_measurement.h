#ifndef ABS_INLINE_MEASUREMENT_H
#define ABS_INLINE_MEASUREMENT_H

#include "chromabase.h"


namespace Chroma { 

  class AbsInlineMeasurement {
  public:

    // Virtual Destructor
    virtual ~AbsInlineMeasurement(void) {}

    // Tell me how often I should measure this beastie
    virtual unsigned long getFrequency(void) const = 0;

    // Do the measurement
    virtual void operator()(const multi1d<LatticeColorMatrix>& u,
			    unsigned long update_no,
			    XMLWriter& xml_out) = 0;
  };

}; // End namespace

#endif

