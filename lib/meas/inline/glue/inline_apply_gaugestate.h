// -*- C++ -*-
// $Id: inline_apply_gaugestate.h,v 3.2 2007-08-27 20:00:22 uid3790 Exp $
/*! \file
 *  \brief Inline gauge state application
 */

#ifndef __inline_apply_gaugestate_h__
#define __inline_apply_gaugestate_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlineGaugeStateEnv 
  {
    extern const std::string name;
    bool registerAll();

    /*! \ingroup inlineglue */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct Param_t
      {
	GroupXML_t    cgs;      /*!< Gauge State */
      } param;

      struct NamedObject_t
      {
	std::string   gauge_id;
	std::string   output_id;
      } named_obj;
    };


    /*! \ingroup inlineglue */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }
}

#endif
