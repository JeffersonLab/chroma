// -*- C++ -*-
// $Id: inline_seqsource_w.h,v 3.3 2006-12-02 18:18:07 edwards Exp $
/*! \file
 * \brief Inline construction of sequential sources
 *
 * Sequential source construction
 */

#ifndef __inline_seqsource_w_h__
#define __inline_seqsource_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineSeqSourceEnv 
  {
    extern const std::string name;
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long      frequency;

      SeqSource_t        param;
      PropSinkSmear_t    sink_header;

      struct NamedObject_t
      {
	std::string            gauge_id;
	multi1d<std::string>   prop_ids;
	std::string            seqsource_id;
      } named_obj;
    };

    //! Compute a sequential source
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }

}

#endif
