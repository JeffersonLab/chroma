// $Id: inline_fuzwilp.cc,v 2.2 2006-02-26 14:17:43 mcneile Exp $
/*! \file
 * \brief Inline fuzzed Wilson loops
 */

#include "chromabase.h"
#include "meas/inline/glue/inline_fuzwilp.h"
#include "meas/glue/fuzwilp.h"
#include "meas/inline/abs_inline_measurement_factory.h"


using namespace QDP;

namespace Chroma 
{ 
  namespace InlineFuzzedWilsonLoopEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineFuzzedWilsonLoop(InlineFuzzedWilsonLoopParams(xml_in, path));
    }
    const std::string name = "FUZZED_WILSON_LOOP";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };

  // Params
  InlineFuzzedWilsonLoopParams::InlineFuzzedWilsonLoopParams() { frequency = 0; }
  InlineFuzzedWilsonLoopParams::InlineFuzzedWilsonLoopParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      read(paramtop, "j_decay", j_decay);
      read(paramtop, "tmax", tmax);
      read(paramtop, "n_smear", n_smear);
      read(paramtop, "sm_fact", sm_fact);
      read(paramtop, "BlkMax",  BlkMax);
      read(paramtop, "BlkAccu", BlkAccu);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }

  void 
  InlineFuzzedWilsonLoop::operator()(const multi1d<LatticeColorMatrix>& u,
				     XMLBufferWriter& gauge_xml,
				     unsigned long update_no,
				     XMLWriter& xml_out) 
  {
    push(xml_out, "APE_Smeared_Wilsonloop");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << "APE_Smeared_Wilsonloop" << endl;

    fuzwilp(u, 
	    params.j_decay, 
            params.tmax,
	    params.n_smear,
	    params.sm_fact,
	    params.BlkAccu,
	    params.BlkMax,
	    xml_out, 
	    "fuzwilp");

    pop(xml_out);
    QDPIO::cout << "FuzzedWilsonLoop ran successfully" << endl;

    END_CODE();
  } 
};
