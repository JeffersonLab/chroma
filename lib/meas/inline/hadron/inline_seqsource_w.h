// -*- C++ -*-
// $Id: inline_seqsource_w.h,v 1.1 2005-04-06 04:34:53 edwards Exp $
/*! \file
 * \brief Inline construction of sequential sources
 *
 * Sequential source construction
 */

#ifndef __inline_seqsource_h__
#define __inline_seqsource_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace InlineSeqSourceEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineSeqSourceParams 
  {
    InlineSeqSourceParams();
    InlineSeqSourceParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long      frequency;

    SeqSource_t        param;
    PropSink_t         sink_header;

    struct Prop_t
    {
      std::string      prop_file;  // The file is expected to be in SciDAC format!
      std::string      seqsource_file;  // The file is expected to be in SciDAC format!
      QDP_volfmt_t     seqsource_volfmt;
    } prop;
  };

  //! Inline measurement of Wilson loops
  class InlineSeqSource : public AbsInlineMeasurement 
  {
  public:
    ~InlineSeqSource() {}
    InlineSeqSource(const InlineSeqSourceParams& p) : params(p) {}
    InlineSeqSource(const InlineSeqSource& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineSeqSourceParams params;
  };

};

#endif
