// -*- C++ -*-
// $Id: inline_random_transf_gauge.h,v 3.1 2007-11-09 21:26:43 edwards Exp $
/*! \file
 *  \brief Do a random gauge transformation on a gauge field
 */

#ifndef __inline_random_transf_gauge_h__
#define __inline_random_transf_gauge_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlineRandomTransfGaugeEnv 
  {
    extern const std::string name;
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineglue */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string     gauge_id;      /*!< input gauge field */
	std::string     rgauge_id;     /*!< output gauge field */
	std::string     gauge_rot_id;  /*!< random gauge transformation fields */
      } named_obj;

    };


    //! Inline random gauge transformation on a gauge field
    /*! \ingroup inlineglue */
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
