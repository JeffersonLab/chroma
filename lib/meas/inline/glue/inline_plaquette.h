// $Id: inline_plaquette.h,v 1.5 2005-04-06 04:34:53 edwards Exp $
/*! \file
 *  \brief Inline plaquette
 */

#ifndef __inline_plaquette_h__
#define __inline_plaquette_h__

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

    unsigned long getFrequency(void) const {return frequency;}

    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    const unsigned long frequency;
  };

};

#endif
