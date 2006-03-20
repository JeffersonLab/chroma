// $Id: inline_fuzwilp.cc,v 2.3 2006-03-20 04:22:02 edwards Exp $
/*! \file
 * \brief Inline fuzzed Wilson loops
 */

#include "chromabase.h"
#include "meas/inline/glue/inline_fuzwilp.h"
#include "meas/glue/fuzwilp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/default_gauge_field.h"


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
  InlineFuzzedWilsonLoopParams::InlineFuzzedWilsonLoopParams()
  { 
    frequency = 0; 
    named_obj.gauge_id = InlineDefaultGaugeField::getId();
  }

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
 
      // Ids
      named_obj.gauge_id = InlineDefaultGaugeField::readGaugeId(paramtop, "NamedObject/gauge_id");
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }

  void 
  InlineFuzzedWilsonLoop::operator()(unsigned long update_no,
				     XMLWriter& xml_out) 
  {
    START_CODE();

    QDP::StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Grab the gauge field
    XMLBufferWriter gauge_xml;
    multi1d<LatticeColorMatrix> u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
    TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

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
 
    snoop.stop();
    QDPIO::cout << InlineFuzzedWilsonLoopEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineFuzzedWilsonLoopEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 
};
