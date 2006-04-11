// $Id: inline_wilslp.cc,v 3.1 2006-04-11 04:18:23 edwards Exp $
/*! \file
 *  \brief Inline Wilson loops
 */

#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/wilslp.h"
#include "meas/inline/io/named_objmap.h"

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



  //! WilsonLoop input
  void read(XMLReader& xml, const string& path, InlineWilsonLoopParams::Param_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "kind", input.kind);
    read(inputtop, "j_decay", input.j_decay);
  }

  //! WilsonLoop output
  void write(XMLWriter& xml, const string& path, const InlineWilsonLoopParams::Param_t& input)
  {
    push(xml, path);

    write(xml, "kind", input.kind);
    write(xml, "j_decay", input.j_decay);

    pop(xml);
  }


  //! WilsonLoop input
  void read(XMLReader& xml, const string& path, InlineWilsonLoopParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! WilsonLoop output
  void write(XMLWriter& xml, const string& path, const InlineWilsonLoopParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  // Param stuff
  InlineWilsonLoopParams::InlineWilsonLoopParams()
  { 
    frequency = 0; 
  }

  InlineWilsonLoopParams::InlineWilsonLoopParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      read(paramtop, "Frequency", frequency);

      // Params
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
    wilslp(u, params.param.j_decay, Nd-1, params.param.kind,
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
