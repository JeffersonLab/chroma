// -*- C++ -*-
// $Id: inline_qio_write_erase_obj.h,v 1.1 2005-09-24 21:14:28 edwards Exp $
/*! \file
 * \brief Inline task to write and delete an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_qio_write_erase_obj_h__
#define __inline_qio_write_erase_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineQIOWriteEraseNamedObjEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineQIOWriteEraseNamedObjParams 
  {
    InlineQIOWriteEraseNamedObjParams();
    InlineQIOWriteEraseNamedObjParams(XMLReader& xml_in, const std::string& path);
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
      QDP_volfmt_t  file_volfmt;
    } file;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineQIOWriteEraseNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineQIOWriteEraseNamedObj() {}
    InlineQIOWriteEraseNamedObj(const InlineQIOWriteEraseNamedObjParams& p) : params(p) {}
    InlineQIOWriteEraseNamedObj(const InlineQIOWriteEraseNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineQIOWriteEraseNamedObjParams params;
  };

};

#endif
