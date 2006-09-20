// -*- C++ -*-
// $Id: inline_list_obj.h,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_list_obj_h__
#define __inline_list_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineListNamedObjEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineListNamedObjParams 
  {
    InlineListNamedObjParams();
    InlineListNamedObjParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;
  };

  //! Inline writing of memory objects
  /*! \ingroup inlineio */
  class InlineListNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineListNamedObj() {}
    InlineListNamedObj(const InlineListNamedObjParams& p) : params(p) {}
    InlineListNamedObj(const InlineListNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineListNamedObjParams params;
  };

};

#endif
