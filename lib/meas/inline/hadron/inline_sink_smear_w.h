// -*- C++ -*-
// $Id: inline_sink_smear_w.h,v 2.1 2005-11-30 04:46:39 edwards Exp $
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
  /*! \ingroup inlinehadron */
  namespace InlineSinkSmearEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineSinkSmearParams 
  {
    InlineSinkSmearParams();
    InlineSinkSmearParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long   frequency;

    PropSinkSmear_t param;

    struct NamedObject_t
    {
      std::string   prop_id;
      std::string   smeared_prop_id;
    } named_obj;
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
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
