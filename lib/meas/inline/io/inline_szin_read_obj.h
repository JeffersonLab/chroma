// -*- C++ -*-
// $Id: inline_szin_read_obj.h,v 2.1 2005-11-01 22:00:01 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_szin_read_obj_h__
#define __inline_szin_read_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineSZINReadNamedObjEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineSZINReadNamedObjParams 
  {
    InlineSZINReadNamedObjParams();
    InlineSZINReadNamedObjParams(XMLReader& xml_in, const std::string& path);
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
  class InlineSZINReadNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineSZINReadNamedObj() {}
    InlineSZINReadNamedObj(const InlineSZINReadNamedObjParams& p) : params(p) {}
    InlineSZINReadNamedObj(const InlineSZINReadNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineSZINReadNamedObjParams params;
  };

};

#endif
