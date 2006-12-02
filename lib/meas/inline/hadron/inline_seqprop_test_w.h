// -*- C++ -*-
// $Id: inline_seqprop_test_w.h,v 3.2 2006-12-02 18:18:07 edwards Exp $
/*! \file
 * \brief Test sequential propagator
 *
 * Sequential source test
 */

#ifndef __inline_seqprop_test_w_h__
#define __inline_seqprop_test_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineSeqPropTestEnv 
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

      PropSourceSmear_t  smear_header;    /*!< Smearing to apply on seqprop at the source */

      //! Propagators
      struct NamedObject_t
      {
	std::string   gauge_id;           /*!< Input Gauge id */
	multi1d<std::string>   sink_ids;  /*!< forward sink smeared propagators needed for 2-pt function */
	std::string   seqprop_id;         /*!< backward propagator */
	int           gamma_insertion;    /*!< second gamma insertion */
      } named_obj;

      std::string xml_file;  /*!< Alternate XML file pattern */
    };

    //! Inline test of sequential propagators
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

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };

  }

}

#endif
