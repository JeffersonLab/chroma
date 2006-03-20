// -*- C++ -*-
// $Id: inline_plaquette.h,v 2.1 2006-03-20 04:22:02 edwards Exp $
/*! \file
 *  \brief Inline plaquette
 */

#ifndef __inline_plaquette_h__
#define __inline_plaquette_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlinePlaquetteEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  /*! \ingroup inlineglue */
  struct InlinePlaquetteParams 
  {
    InlinePlaquetteParams();
    InlinePlaquetteParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;

    struct NamedObject_t
    {
      std::string   gauge_id;
    } named_obj;
  };


  /*! \ingroup inlineglue */
  class InlinePlaquette : public AbsInlineMeasurement 
  {
  public:
    ~InlinePlaquette() {}
    InlinePlaquette(const InlinePlaquetteParams& p) : params(p) {}
    InlinePlaquette(const InlinePlaquette& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    void operator()(unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlinePlaquetteParams params;
  };

};

#endif
