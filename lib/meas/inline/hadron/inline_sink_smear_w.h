// -*- C++ -*-
// $Id: inline_sink_smear_w.h,v 3.3 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline sink_smear propagators
 *
 * Sink smear propagators
 */

#ifndef __inline_sink_smear_w_h__
#define __inline_sink_smear_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineSinkSmearEnv 
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

      unsigned long   frequency;

      PropSinkSmear_t param;

      struct NamedObject_t
      {
	std::string   gauge_id;
	std::string   prop_id;
	std::string   smeared_prop_id;
      } named_obj;
    };


    //! Inline task for sinking smearing propagators
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
