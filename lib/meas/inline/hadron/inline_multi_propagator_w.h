// -*- C++ -*-
// $Id: inline_multi_propagator_w.h,v 1.1 2005-07-15 11:06:11 bjoo Exp $
/*! \file
 *  \brief Inline construction of multi_propagator -- overlap only
 *  
 * 
 *
 * Propagator calculations
 */

#ifndef __inline_multi_propagator_h__
#define __inline_multi_propagator_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMultiPropagatorEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */ 
  struct InlineMultiPropagatorParams 
  {
    InlineMultiPropagatorParams();
    InlineMultiPropagatorParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long     frequency;

    ChromaMultiProp_t      param;
    std::string            stateInfo;

    struct Prop_t
    {
      std::string     source_file;
      std::string     prop_file;
      QDP_volfmt_t    prop_volfmt;
    } prop;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineMultiPropagator : public AbsInlineMeasurement 
  {
  public:
    ~InlineMultiPropagator() {}
    InlineMultiPropagator(const InlineMultiPropagatorParams& p) : params(p) {}
    InlineMultiPropagator(const InlineMultiPropagator& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const multi1d<LatticeColorMatrix>& u,
	      XMLBufferWriter& gauge_xml,
	      const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineMultiPropagatorParams params;
  };

};

#endif
