// $Id: inline_wilslp.cc,v 1.2 2005-04-06 04:34:53 edwards Exp $
/*! \file
 *  \brief Inline Wilson loops
 */

#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/wilslp.h"

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
  InlineWilsonLoopParams::InlineWilsonLoopParams() { frequency = 0; }

  InlineWilsonLoopParams::InlineWilsonLoopParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      read(paramtop, "Frequency", frequency);
      read(paramtop, "kind", kind);
      read(paramtop, "j_decay", j_decay);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void 
  InlineWilsonLoop::operator()(const multi1d<LatticeColorMatrix>& u,
			       XMLBufferWriter& gauge_xml,
			       unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    push(xml_out, "WilsonLoop");
    write(xml_out, "update_no", update_no);

    Double w_plaq, s_plaq, t_plaq, link;
    wilslp(u, params.j_decay, Nd-1, params.kind,
	   xml_out, "wilslp");
    
    pop(xml_out); // pop("WilsonLoop");
    
  } 

};
