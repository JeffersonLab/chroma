// -*- C++ -*-
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_qio_write_obj_h__
#define __inline_qio_write_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineQIOWriteNamedObjEnv 
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
	std::string   object_type;
      } named_obj;

      struct File_t
      {
	std::string   file_name;
	QDP_volfmt_t  file_volfmt;
	bool          parallel_io;
      } file;
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
