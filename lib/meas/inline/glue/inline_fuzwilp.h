// -*- C++ -*-
// $Id: inline_fuzwilp.h,v 3.2 2006-09-20 20:28:01 edwards Exp $
/*! \file
 * \brief Inline fuzzed Wilson loops
 */

#ifndef __inline_fuzwilp_h__
#define __inline_fuzwilp_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlineFuzzedWilsonLoopEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Parameter structure
  /*! \ingroup inlineglue */
  struct InlineFuzzedWilsonLoopParams 
  {
    InlineFuzzedWilsonLoopParams();
    InlineFuzzedWilsonLoopParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int j_decay;
      int tmax;
      int n_smear;
      int BlkMax;
      Real sm_fact;
      Real BlkAccu; 
    } param;

    struct NamedObject_t
    {
      std::string   gauge_id;
    } named_obj;
  };

  //! Inline measurement of fuzzed Wilson loops
  /*! \ingroup inlineglue */
  class InlineFuzzedWilsonLoop : public AbsInlineMeasurement 
  {
  public:
    ~InlineFuzzedWilsonLoop() {}
    InlineFuzzedWilsonLoop(const InlineFuzzedWilsonLoopParams& p) : params(p) {}
    InlineFuzzedWilsonLoop(const InlineFuzzedWilsonLoop& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineFuzzedWilsonLoopParams params;
  };

};

#endif
