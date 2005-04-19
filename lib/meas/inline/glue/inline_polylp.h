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
    InlinePolyakovLoopParams() { frequency = 0; }

    InlinePolyakovLoopParams(XMLReader& xml_in, const std::string& path) 
    {
      try {
	XMLReader paramtop(xml_in, path);
	read(paramtop, "./Frequency", frequency);
      }
      catch(const std::string& e) { 
	QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }

    unsigned long frequency;
  };


  /*! \ingroup inlineglue */
  class InlinePolyakovLoop : public AbsInlineMeasurement 
  {
  public:
    ~InlinePolyakovLoop() {}
    InlinePolyakovLoop(const InlinePolyakovLoopParams& p_) : p(p_) {}
    InlinePolyakovLoop(const InlinePolyakovLoop& p_) : p(p_.p) {}

    unsigned long getFrequency(void) const {return p.frequency;}

    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    const InlinePolyakovLoopParams p;
  };

};

#endif
