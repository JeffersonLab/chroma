// -*- C++ -*-
// $Id: inline_propagator_s.h,v 3.1 2006-11-17 02:17:31 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#ifndef __inline_propagator_s_h__
#define __inline_propagator_s_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineStaggeredPropagatorEnv 
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

      unsigned long     frequency;

      ChromaProp_t      param;

      struct NamedObject_t
      {
	std::string     gauge_id;
	std::string     source_id;
	std::string     prop_id;
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };

    //! Inline computation of propagators
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
