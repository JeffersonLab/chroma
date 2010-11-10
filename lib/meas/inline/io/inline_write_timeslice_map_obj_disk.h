// -*- C++ -*-
/*! \file
 * \brief Inline task to time-sliced map object
 */

#ifndef __inline_write_timeslice_map_obj_disk_h__
#define __inline_write_timeslice_map_obj_disk_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineWriteTimeSliceMapObjDiskEnv 
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
	std::string   object_type;         /*!< Input object type */
	std::string   input_id;            /*!< Input object id */
	std::string   output_file;         /*!< Output map-object-disk */
      };

      NamedObject_t   named_obj;
    };


    //! Inline task to time-sliced map object
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
