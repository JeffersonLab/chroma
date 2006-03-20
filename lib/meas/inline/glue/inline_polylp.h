// -*- C++ -*-

#ifndef INLINE_POLYLOOP_LOOP_H
#define INLINE_POLYAKOV_LOOP_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"



namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlinePolyakovLoopEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  /*! \ingroup inlineglue */
  struct InlinePolyakovLoopParams 
  {
    InlinePolyakovLoopParams();
    InlinePolyakovLoopParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;

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
