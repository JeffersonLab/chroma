// -*- C++ -*-
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_qio_read_obj_h__
#define __inline_qio_read_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineQIOReadNamedObjEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string   object_id;
	std::string   object_type;
      };

      struct File_t
      {
	std::string   file_name;
	bool parallel_io;
      };

      File_t file;
      NamedObject_t named_obj;
      GroupXML_t    named_obj_xml;  /*!< Holds standard named objects */
    };

    //! Inline reading of qio objects
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

  } // namespace InlineQIOReadNamedObjEnv 

} // namespace Chroma

#endif
