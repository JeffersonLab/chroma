#ifndef INLINE_PLAQUETTE_H
#define INLINE_PLAQUETTE_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"



namespace Chroma { 
  namespace InlinePlaquetteEnv {
    extern const std::string name;
    extern const bool registered;
  }

  struct InlinePlaquetteParams {
    InlinePlaquetteParams() { frequency = 0; }

    InlinePlaquetteParams(XMLReader& xml_in, const std::string& path) 
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

  class InlinePlaquette : public AbsInlineMeasurement {
  public:
    ~InlinePlaquette() {}

    InlinePlaquette(const InlinePlaquetteParams& p) : frequency(p.frequency) {}

    InlinePlaquette(const InlinePlaquette& p) : frequency(p.frequency) {}

    const unsigned long getFrequency(void) const {
      return frequency;
    }


    void operator()(const multi1d<LatticeColorMatrix>& u,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    const unsigned long frequency;
  };

};

#endif
