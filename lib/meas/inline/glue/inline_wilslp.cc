// $Id: inline_wilslp.cc,v 3.3 2006-07-10 19:04:47 edwards Exp $
/*! \file
 *  \brief Inline Wilson loops
 */

#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/wilslp.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/gauge/gaugebcs/gaugebc_factory.h"
#include "actions/gauge/gaugebcs/gaugebc_aggregate.h"

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

    bool registerAll()
    {
      bool foo = true;
      foo &= GaugeTypeGaugeBCEnv::registered;
      foo &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
      return foo;
    }

    const bool registered = registerAll();
  };



  //! WilsonLoop input
  void read(XMLReader& xml, const string& path, InlineWilsonLoopParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 2:
      param.gaugebc = readXMLGroup(paramtop, "GaugeBC", "Name");
      break;

    default:
      QDPIO::cerr << "InlineWilsonLoopParams::Param_t: " << version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "kind", param.kind);
    read(paramtop, "j_decay", param.j_decay);
    read(paramtop, "t_dir", param.t_dir);
  }

  //! WilsonLoop output
  void write(XMLWriter& xml, const string& path, const InlineWilsonLoopParams::Param_t& param)
  {
    push(xml, path);

    int version = 2;
    write(xml, "version", version);
    write(xml, "kind", param.kind);
    write(xml, "j_decay", param.j_decay);
    write(xml, "t_dir", param.t_dir);
    xml << param.gaugebc.xml;

    pop(xml);
  }


  //! WilsonLoop input
  void read(XMLReader& xml, const string& path, InlineWilsonLoopParams::NamedObject_t& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "gauge_id", param.gauge_id);
  }

  //! WilsonLoop output
  void write(XMLWriter& xml, const string& path, const InlineWilsonLoopParams::NamedObject_t& param)
  {
    push(xml, path);

    write(xml, "gauge_id", param.gauge_id);

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

    try
    {
      // Grab the gauge field
      multi1d<LatticeColorMatrix> u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      // Set the gaugebc and modify the fields
      {
	std::istringstream  xml_s(params.param.gaugebc.xml);
	XMLReader  gaugetop(xml_s);

	Handle<GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  gbc(TheGaugeBCFactory::Instance().createObject(params.param.gaugebc.id,
							 gaugetop,
							 params.param.gaugebc.path));
	gbc->zero(u);
      }
    
      push(xml_out, "WilsonLoop");
      write(xml_out, "update_no", update_no);
      write(xml_out, "decay_dir", params.param.j_decay);
      write(xml_out, "t_dir", params.param.t_dir);

      Double w_plaq, s_plaq, t_plaq, link;
      wilslp(u, params.param.j_decay, params.param.t_dir, params.param.kind,
	     xml_out, "wilslp");
    
      pop(xml_out); // pop("WilsonLoop");
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineWilsonLoopEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineWilsonLoopEnv::name << ": caught error: " << e << endl;
      QDP_abort(1);
    }
 
    snoop.stop();
    QDPIO::cout << InlineWilsonLoopEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineWilsonLoopEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
