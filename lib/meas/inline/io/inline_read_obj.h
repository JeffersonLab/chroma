// -*- C++ -*-
// $Id: inline_read_obj.h,v 1.1 2005-09-23 03:43:09 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_read_obj_h__
#define __inline_read_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineReadNamedObjEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineReadNamedObjParams 
  {
    InlineReadNamedObjParams();
    InlineReadNamedObjParams(XMLReader& xml_in, const std::string& path);
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
  class InlineReadNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineReadNamedObj() {}
    InlineReadNamedObj(const InlineReadNamedObjParams& p) : params(p) {}
    InlineReadNamedObj(const InlineReadNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineReadNamedObjParams params;
  };

};

#endif
