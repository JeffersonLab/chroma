// -*- C++ -*-
// $Id: inline_propagator_w.h,v 1.1 2005-04-06 04:34:53 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#ifndef __inline_propagator_h__
#define __inline_propagator_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace InlinePropagatorEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlinePropagatorParams 
  {
    InlinePropagatorParams();
    InlinePropagatorParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    ChromaProp_t      param;
    std::string       stateInfo;

    struct Prop_t
    {
      std::string     source_file;
      std::string     prop_file;
      QDP_volfmt_t    prop_volfmt;
    } prop;
  };

  //! Inline measurement of Wilson loops
  class InlinePropagator : public AbsInlineMeasurement 
  {
  public:
    ~InlinePropagator() {}
    InlinePropagator(const InlinePropagatorParams& p) : params(p) {}
    InlinePropagator(const InlinePropagator& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlinePropagatorParams params;
  };

};

#endif
