// -*- C++ -*-
// $Id: inline_snarf_obj.h,v 1.1 2005-09-24 21:14:28 edwards Exp $
/*! \file
 * \brief Inline task to copy gauge arg to a named object
 *
 * Named object writing
 */

#ifndef __inline_snarf_obj_h__
#define __inline_snarf_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineSnarfNamedObjEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlineio */
  struct InlineSnarfNamedObjParams 
  {
    InlineSnarfNamedObjParams();
    InlineSnarfNamedObjParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct NamedObject_t
    {
      std::string   object_id;
    } named_obj;
  };

  //! Inline snarfing of memory objects
  /*! \ingroup inlineio */
  class InlineSnarfNamedObj : public AbsInlineMeasurement 
  {
  public:
    ~InlineSnarfNamedObj() {}
    InlineSnarfNamedObj(const InlineSnarfNamedObjParams& p) : params(p) {}
    InlineSnarfNamedObj(const InlineSnarfNamedObj& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the writing
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineSnarfNamedObjParams params;
  };

};

#endif
