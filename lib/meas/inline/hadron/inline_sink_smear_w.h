// -*- C++ -*-
// $Id: inline_sink_smear_w.h,v 1.1 2005-04-06 04:34:54 edwards Exp $
/*! \file
 * \brief Inline sink_smear propagators
 *
 * Sink smear propagators
 */

#ifndef __inline_sink_smear_h__
#define __inline_sink_smear_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace InlineSinkSmearEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineSinkSmearParams 
  {
    InlineSinkSmearParams();
    InlineSinkSmearParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long   frequency;

    PropSink_t      param;

    struct Prop_t
    {
      std::string   prop_file;

      QDP_volfmt_t  smeared_prop_volfmt; // Volume format is  SINGLEFILE or MULTIFILE
      std::string   smeared_prop_file;
    } prop;
  };


  //! Inline measurement of Wilson loops
  class InlineSinkSmear : public AbsInlineMeasurement 
  {
  public:
    ~InlineSinkSmear() {}
    InlineSinkSmear(const InlineSinkSmearParams& p) : params(p) {}
    InlineSinkSmear(const InlineSinkSmear& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineSinkSmearParams params;
  };

};

#endif
