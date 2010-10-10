// -*- C++ -*-
/*! \file
 * \brief Inline task to read an object into a named buffer
 */

#ifndef __inline_read_map_obj_disk_h__
#define __inline_read_map_obj_disk_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineReadMapObjDiskEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params(XMLReader& xml_in, const std::string& path);

      unsigned int frequency;

      struct File {
	std::string   file_name;
      } file;

      struct NamedObject_t {
	std::string   object_type; 
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
