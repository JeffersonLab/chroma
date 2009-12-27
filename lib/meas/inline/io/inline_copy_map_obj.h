// -*- C++ -*-
/*! \file
 * \brief Inline task to copy map objects
 */

#ifndef __inline_copy_map_obj_h__
#define __inline_copy_map_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineCopyMapObjEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned int frequency;

      struct NamedObject_t {
	std::string   object_type;
	std::string   input_id;
	std::string   output_id;
	GroupXML_t    output_obj;  /*!< Output object map */
      };

      NamedObject_t   named_obj;
    };


    //! Inline copying of map objects
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
