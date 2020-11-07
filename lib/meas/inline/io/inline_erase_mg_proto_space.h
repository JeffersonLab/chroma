// -*- C++ -*-
/*! \file
 * \brief Inline task to erase a named mg space
 *
 * Named object writing
 */

#ifndef __inline_erase_mg_proto_space_h__
#define __inline_erase_mg_proto_space_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"


namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineEraseMGProtoSpaceEnv
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string   object_id;
      } named_obj;
    };

    //! Inline writing of memory objects
    /*! \ingroup inlineio */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }

}

#endif
