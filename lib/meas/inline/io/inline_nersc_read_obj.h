// -*- C++ -*-
// $Id: inline_nersc_read_obj.h,v 3.3 2008-06-18 21:38:28 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_nersc_read_obj_h__
#define __inline_nersc_read_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
//#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineNERSCReadNamedObjEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineNERSCReadNamedObjParams 
  {
    InlineNERSCReadNamedObjParams();
    InlineNERSCReadNamedObjParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct NamedObject_t
    {
      std::string   object_id;
    } named_obj;

    struct File_t
    {
      std::string   file_name;
    } file;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineNERSCReadNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineNERSCReadNamedObj() {}
    InlineNERSCReadNamedObj(const InlineNERSCReadNamedObjParams& p) : params(p) {}
    InlineNERSCReadNamedObj(const InlineNERSCReadNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineNERSCReadNamedObjParams params;
  };

};

#endif
