// -*- C++ -*-
// $Id: inline_wilslp.h,v 1.1 2005-02-10 15:50:47 edwards Exp $
/*! \file
 *  \brief Inline Wilson loops
 */

#ifndef __inline_wilslp_h__
#define __inline_wilslp_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  namespace InlineWilsonLoopEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineWilsonLoopParams 
  {
    InlineWilsonLoopParams();
    InlineWilsonLoopParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;
    int kind;
    int j_decay;
  };

  //! Inline measurement of Wilson loops
  class InlineWilsonLoop : public AbsInlineMeasurement 
  {
  public:
    ~InlineWilsonLoop() {}
    InlineWilsonLoop(const InlineWilsonLoopParams& p) : params(p) {}
    InlineWilsonLoop(const InlineWilsonLoop& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineWilsonLoopParams params;
  };

};

#endif
