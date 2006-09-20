// $Id: inline_fuzwilp.cc,v 3.2 2006-09-20 20:28:01 edwards Exp $
/*! \file
 * \brief Inline fuzzed Wilson loops
 */

#include "chromabase.h"
#include "meas/inline/glue/inline_fuzwilp.h"
#include "meas/glue/fuzwilp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"


using namespace QDP;

namespace Chroma 
{ 
  namespace InlineFuzzedWilsonLoopEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineFuzzedWilsonLoop(InlineFuzzedWilsonLoopParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "FUZZED_WILSON_LOOP";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }


  //! FuzzedWilsonLoop input
  void read(XMLReader& xml, const string& path, InlineFuzzedWilsonLoopParams::Param_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "j_decay", input.j_decay);
    read(inputtop, "tmax", input.tmax);
    read(inputtop, "n_smear", input.n_smear);
    read(inputtop, "sm_fact", input.sm_fact);
    read(inputtop, "BlkMax",  input.BlkMax);
    read(inputtop, "BlkAccu", input.BlkAccu);
   }

  //! FuzzedWilsonLoop output
  void write(XMLWriter& xml, const string& path, const InlineFuzzedWilsonLoopParams::Param_t& input)
  {
    push(xml, path);

    write(xml, "j_decay", input.j_decay);
    write(xml, "tmax", input.tmax);
    write(xml, "n_smear", input.n_smear);
    write(xml, "sm_fact", input.sm_fact);
    write(xml, "BlkMax",  input.BlkMax);
    write(xml, "BlkAccu", input.BlkAccu);

    pop(xml);
  }


  //! FuzzedWilsonLoop input
  void read(XMLReader& xml, const string& path, InlineFuzzedWilsonLoopParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! FuzzedWilsonLoop output
  void write(XMLWriter& xml, const string& path, const InlineFuzzedWilsonLoopParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  // Params
  InlineFuzzedWilsonLoopParams::InlineFuzzedWilsonLoopParams()
  { 
    frequency = 0; 
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

      // params
      read(paramtop, "Param", param);

      // Ids
      read(paramtop, "NamedObject", named_obj);
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
	    params.param.j_decay, 
            params.param.tmax,
	    params.param.n_smear,
	    params.param.sm_fact,
	    params.param.BlkAccu,
	    params.param.BlkMax,
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
