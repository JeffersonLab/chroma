// -*- C++ -*-
// $Id: inline_xml_write_obj.h,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_xml_write_obj_h__
#define __inline_xml_write_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineXMLWriteNamedObjEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineXMLWriteNamedObjParams 
  {
    InlineXMLWriteNamedObjParams();
    InlineXMLWriteNamedObjParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct NamedObject_t
    {
      std::string   object_id;
      std::string   object_type;
    } named_obj;

    struct File_t
    {
      std::string   file_name;
    } file;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineXMLWriteNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineXMLWriteNamedObj() {}
    InlineXMLWriteNamedObj(const InlineXMLWriteNamedObjParams& p) : params(p) {}
    InlineXMLWriteNamedObj(const InlineXMLWriteNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineXMLWriteNamedObjParams params;
  };

};

#endif
