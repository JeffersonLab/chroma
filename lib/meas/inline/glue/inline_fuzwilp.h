// -*- C++ -*-
// $Id: inline_fuzwilp.h,v 2.0 2005-09-25 21:04:37 edwards Exp $
/*! \file
 * \brief Inline spectrum calculations
 *
 * Spectrum calculations
 */

#ifndef __inline_fuzwilp_h__
#define __inline_fuzwilp_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  namespace InlineFuzzedWilsonLoopEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  struct InlineFuzzedWilsonLoopParams 
  {
    InlineFuzzedWilsonLoopParams();
    InlineFuzzedWilsonLoopParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;
    int j_decay;
    int n_smear;
    int BlkMax;
    Real sm_fact;
    Real BlkAccu;
  };

  class InlineFuzzedWilsonLoop : public AbsInlineMeasurement 
  {
  public:
    ~InlineFuzzedWilsonLoop() {}
    InlineFuzzedWilsonLoop(const InlineFuzzedWilsonLoopParams& p) : params(p) {}
    InlineFuzzedWilsonLoop(const InlineFuzzedWilsonLoop& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineFuzzedWilsonLoopParams params;
  };

};

#endif
