#ifndef INLINE_POLYLOOP_LOOP_H
#define INLINE_POLYAKOV_LOOP_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"


using namespace QDP;
using namespace Chroma;

namespace Chroma { 
  namespace InlinePolyakovLoopEnv {
    extern const std::string name;
    extern const bool registered;
  }

  struct InlinePolyakovLoopParams {
    InlinePolyakovLoopParams() { frequency = 0; }

    InlinePolyakovLoopParams(XMLReader& xml_in, const std::string& path) 
    {
      try {
	XMLReader paramtop(xml_in, path);
	read(paramtop, "./Frequency", frequency);
      }
      catch(const std::string& e) { 
	QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }

    unsigned long frequency;
  };

  class InlinePolyakovLoop : public AbsInlineMeasurement {
  public:
    ~InlinePolyakovLoop() {}

    InlinePolyakovLoop(const InlinePolyakovLoopParams& p_) : p(p_) {}

    InlinePolyakovLoop(const InlinePolyakovLoop& p_) : p(p_.p) {}

    const unsigned long getFrequency(void) const {
      return p.frequency;
    }


    void operator()(const multi1d<LatticeColorMatrix>& u,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    const InlinePolyakovLoopParams p;
  };

};

#endif
