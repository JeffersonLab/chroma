// -*- C++ -*-
// $Id: inline_make_source_w.h,v 1.2 2005-04-19 20:05:22 edwards Exp $
/*! \file
 * \brief Inline construction of make_source
 *
 * Construct source for propagator calculations
 */

#ifndef __inline_make_source_h__
#define __inline_make_source_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMakeSourceEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMakeSourceParams 
  {
    InlineMakeSourceParams();
    InlineMakeSourceParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    PropSource_t  param;

    struct Prop_t
    {
      std::string     source_file;
      QDP_volfmt_t    source_volfmt;
    } prop;
  };

  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineMakeSource : public AbsInlineMeasurement 
  {
  public:
    ~InlineMakeSource() {}
    InlineMakeSource(const InlineMakeSourceParams& p) : params(p) {}
    InlineMakeSource(const InlineMakeSource& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineMakeSourceParams params;
  };

};

#endif
