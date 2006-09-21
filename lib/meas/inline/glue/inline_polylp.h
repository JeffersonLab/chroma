// -*- C++ -*-
// $Id: inline_polylp.h,v 3.2 2006-09-21 18:43:27 edwards Exp $
/*! \file
 *  \brief Inline polyakov loop
 */

#ifndef INLINE_POLYLOOP_LOOP_H
#define INLINE_POLYAKOV_LOOP_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlinePolyakovLoopEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! \ingroup inlineglue */
  struct InlinePolyakovLoopParams 
  {
    InlinePolyakovLoopParams();
    InlinePolyakovLoopParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      GroupXML_t    cgs;      /*!< Gauge State */
    } param;

    struct NamedObject_t
    {
      std::string   gauge_id;
    } named_obj;
  };


  /*! \ingroup inlineglue */
  class InlinePolyakovLoop : public AbsInlineMeasurement 
  {
  public:
    ~InlinePolyakovLoop() {}
    InlinePolyakovLoop(const InlinePolyakovLoopParams& p) : params(p) {}
    InlinePolyakovLoop(const InlinePolyakovLoop& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    void operator()(unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlinePolyakovLoopParams params;
  };

};

#endif
