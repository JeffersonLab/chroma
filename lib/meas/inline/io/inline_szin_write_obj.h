// -*- C++ -*-
// $Id: inline_szin_write_obj.h,v 3.0 2006-04-03 04:59:03 edwards Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_szin_write_obj_h__
#define __inline_szin_write_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineSZINWriteNamedObjEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineSZINWriteNamedObjParams 
  {
    InlineSZINWriteNamedObjParams();
    InlineSZINWriteNamedObjParams(XMLReader& xml_in, const std::string& path);
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

      /*
       *  Now some various rules for truncating the configuration
       */
      bool          trunc;	// Whether to truncate the output
      int           j_decay;    // Direction of time
      int           t_start;	// Starting time slice
      int           t_end;	// Ending time slice
    } file;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineSZINWriteNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineSZINWriteNamedObj() {}
    InlineSZINWriteNamedObj(const InlineSZINWriteNamedObjParams& p) : params(p) {}
    InlineSZINWriteNamedObj(const InlineSZINWriteNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineSZINWriteNamedObjParams params;
  };

};

#endif
