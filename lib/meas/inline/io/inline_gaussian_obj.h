// -*- C++ -*-
// $Id: inline_gaussian_obj.h,v 3.0 2006-04-03 04:59:03 edwards Exp $
/*! \file
 * \brief Inline task to gaussian init a named object
 *
 * Named object initialization
 */

#ifndef __inline_gaussian_obj_h__
#define __inline_gaussian_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineGaussianInitNamedObjEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineGaussianInitNamedObjParams 
  {
    InlineGaussianInitNamedObjParams();
    InlineGaussianInitNamedObjParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct NamedObject_t
    {
      std::string   object_id;
      std::string   object_type;
    } named_obj;
  };

  //! Inline gaussianing of memory objects
  /*! \ingroup inlineio */
  class InlineGaussianInitNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineGaussianInitNamedObj() {}
    InlineGaussianInitNamedObj(const InlineGaussianInitNamedObjParams& p) : params(p) {}
    InlineGaussianInitNamedObj(const InlineGaussianInitNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineGaussianInitNamedObjParams params;
  };

};

#endif
