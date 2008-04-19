// $Id: inline_wilslp.cc,v 3.12 2008-04-19 03:13:32 edwards Exp $
/*! \file
 *  \brief Inline Wilson loops
 */

#include "meas/inline/glue/inline_wilslp.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/glue/wilslp.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/make_xml_file.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma 
{ 
  namespace InlineWilsonLoopEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineWilsonLoop(InlineWilsonLoopParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "WILSLP";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateGaugeStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  //! WilsonLoop input
  void read(XMLReader& xml, const string& path, InlineWilsonLoopParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    param.cgs = CreateGaugeStateEnv::nullXMLGroup();

    switch (version) 
    {
    case 3:
      if (paramtop.count("GaugeState") != 0)
	param.cgs = readXMLGroup(paramtop, "GaugeState", "Name");
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

    int version = 3;
    write(xml, "version", version);
    write(xml, "kind", param.kind);
    write(xml, "j_decay", param.j_decay);
    write(xml, "t_dir", param.t_dir);
    xml << param.cgs.xml;

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

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  // Function call
  void 
  InlineWilsonLoop::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "WilsonLoop");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
    else
    {
      func(update_no, xml_out);
    }
  }


  // Real work done here
  void 
  InlineWilsonLoop::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    QDPIO::cout << InlineWilsonLoopEnv::name << ": Wilson-loop measurements" << endl;

    QDP::StopWatch snoop;
    snoop.reset();
    snoop.start();

    push(xml_out, "WilsonLoop");
    write(xml_out, "update_no", update_no);
    write(xml_out, "decay_dir", params.param.j_decay);
    write(xml_out, "t_dir", params.param.t_dir);

    try
    {
      // Grab the gauge field
      multi1d<LatticeColorMatrix> u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      // Calculate some gauge invariant observables
      MesPlq(xml_out, "Observables", u);

      // Set the construct state and modify the fields
      {
	std::istringstream  xml_s(params.param.cgs.xml);
	XMLReader  gaugetop(xml_s);

	Handle<CreateGaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  cgs(TheCreateGaugeStateFactory::Instance().createObject(params.param.cgs.id,
								  gaugetop,
								  params.param.cgs.path));

	Handle<GaugeState< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  state((*cgs)(u));

	// Pull the u fields back out from the state since they might have been
	// munged with gaugeBC's
	u = state->getLinks();
      }
    
      // Again calculate some gauge invariant observables
      MesPlq(xml_out, "Link_observables", u);

      // Compute the Wilson loops
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
