// -*- C++ -*-
/*! \file
 *  \brief Inline plaquette density
 */

#ifndef INLINE_PLAQ_DENSITY_H
#define INLINE_PLAQ_DENSITY_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlinePlaqDenEnv 
  {
    bool registerAll();

    /*! \ingroup inlineglue */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long frequency;

      struct Param_t
      {
      } param;

      struct NamedObject_t
      {
	std::string   gauge_id;
	std::string   out_file;
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

    protected:
      void func(const unsigned long update_no, 
		XMLWriter& xml_out);

    private:
      Params params;
    };

  }

}

#endif
