// $Id: inline_wilslp.cc,v 2.1 2006-03-20 04:22:02 edwards Exp $
/*! \file
 *  \brief Inline Wilson loops
 */

#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/wilslp.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/default_gauge_field.h"

namespace Chroma 
{ 
  namespace InlineWilsonLoopEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineWilsonLoop(InlineWilsonLoopParams(xml_in, path));
    }

    const std::string name = "WILSLP";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };



  // Param stuff
  InlineWilsonLoopParams::InlineWilsonLoopParams()
  { 
    frequency = 0; 
    named_obj.gauge_id = InlineDefaultGaugeField::getId();
  }

  InlineWilsonLoopParams::InlineWilsonLoopParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      read(paramtop, "Frequency", frequency);
      read(paramtop, "kind", kind);
      read(paramtop, "j_decay", j_decay);

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
  InlineWilsonLoop::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    START_CODE();

    QDP::StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Grab the gauge field
    multi1d<LatticeColorMatrix> u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "WilsonLoop");
    write(xml_out, "update_no", update_no);

    Double w_plaq, s_plaq, t_plaq, link;
    wilslp(u, params.j_decay, Nd-1, params.kind,
	   xml_out, "wilslp");
    
    pop(xml_out); // pop("WilsonLoop");

 
    snoop.stop();
    QDPIO::cout << InlineWilsonLoopEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineWilsonLoopEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
