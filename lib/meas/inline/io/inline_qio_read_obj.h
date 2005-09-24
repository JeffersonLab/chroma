// -*- C++ -*-
// $Id: inline_qio_read_obj.h,v 1.1 2005-09-24 21:14:28 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_qio_read_obj_h__
#define __inline_qio_read_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineQIOReadNamedObjEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineQIOReadNamedObjParams 
  {
    InlineQIOReadNamedObjParams();
    InlineQIOReadNamedObjParams(XMLReader& xml_in, const std::string& path);
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
  class InlineQIOReadNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineQIOReadNamedObj() {}
    InlineQIOReadNamedObj(const InlineQIOReadNamedObjParams& p) : params(p) {}
    InlineQIOReadNamedObj(const InlineQIOReadNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQIOReadNamedObjParams params;
  };

};

#endif
