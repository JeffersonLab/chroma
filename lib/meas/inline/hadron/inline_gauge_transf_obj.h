// -*- C++ -*-
// $Id: inline_gauge_transf_obj.h,v 3.1 2007-11-09 21:26:56 edwards Exp $
/*! \file
 * \brief Inline task gauge transform some fermion object
 *
 * Gauge transform a named object
 */

#ifndef __inline_gauge_transf_obj_h__
#define __inline_gauge_transf_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineGaugeTransfNamedObjEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string   gauge_rot_id;  /*!< random gauge transformation fields */
	std::string   input_id;      /*!< input object */
	std::string   output_id;     /*!< output object */
	std::string   object_type;   /*!< type of the object, like LatticePropagator, etc. */
      } named_obj;
    };

    //! Gauge transform a named object
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // InlineGaugeTransfNamedObjEnv

} // namespace Chroma

#endif
