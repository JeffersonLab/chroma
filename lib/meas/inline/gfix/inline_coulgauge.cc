// $Id: inline_coulgauge.cc,v 3.3 2007-11-09 21:18:09 edwards Exp $
/*! \file
 *  \brief Inline coulomb (and landau) gauge fixing loops
 */

#include "meas/inline/gfix/inline_coulgauge.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/gfix/coulgauge.h"
#include "meas/glue/mesplq.h"
#include "util/info/proginfo.h"
#include "util/gauge/unit_check.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  //! Parameters for running code
  void read(XMLReader& xml, const string& path, InlineCoulGaugeEnv::Params::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      break;

    default :

      QDPIO::cerr << "Input version " << version << " unsupported." << endl;
      QDP_abort(1);
    }
    
    read(paramtop, "j_decay", param.j_decay);
    read(paramtop, "GFAccu", param.GFAccu);
    read(paramtop, "GFMax", param.GFMax);
    read(paramtop, "OrDo", param.OrDo);
    read(paramtop, "OrPara", param.OrPara);
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const InlineCoulGaugeEnv::Params::Param_t& param)
  {
    push(xml, path);
    
    int version = 1;
    write(xml, "version", version);
    write(xml, "GFAccu", param.GFAccu);
    write(xml, "GFMax", param.GFMax);
    write(xml, "OrDo", param.OrDo);
    write(xml, "OrPara", param.OrPara);
    write(xml, "j_decay", param.j_decay);

    pop(xml);
  }

  // Reader for out gauge file
  void read(XMLReader& xml, const string& path, InlineCoulGaugeEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "gfix_id", input.gfix_id);
    read(inputtop, "gauge_rot_id", input.gauge_rot_id);
  }

  // Reader for out gauge file
  void write(XMLWriter& xml, const string& path, const InlineCoulGaugeEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "gfix_id", input.gfix_id);
    write(xml, "gauge_rot_id", input.gauge_rot_id);

    pop(xml);
  }



  namespace InlineCoulGaugeEnv 
  { 
    //! Anonymous namespace
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "COULOMB_GAUGEFIX";

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


    // Param stuff
    Params::Params()
    { 
      frequency = 0; 
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;
      
	// Read program parameters
	read(paramtop, "Param", param);

	// Read in the gfix outfile
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << InlineCoulGaugeEnv::name << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    // Write params
    void
    Params::write(XMLWriter& xml, const std::string& path) 
    {
      push(xml, path);
      
      Chroma::write(xml, "Param", param);
      Chroma::write(xml, "NamedObject", named_obj);

      pop(xml);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      QDP::StopWatch snoop;
      snoop.reset();
      snoop.start();

      XMLBufferWriter gauge_xml;
      multi1d<LatticeColorMatrix> u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

      push(xml_out, "CoulGauge");
      write(xml_out, "update_no", update_no);
    
      QDPIO::cout << InlineCoulGaugeEnv::name << ": coulomb gauge fix" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.write(xml_out, "Input");

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables
      MesPlq(xml_out, "Observables", u);

      // Now coulomb (landau) gauge fix
      multi1d<LatticeColorMatrix> u_gfix(Nd);
      u_gfix = u;

      LatticeColorMatrix g;  // the gauge rotation fields

      int n_gf;
      coulGauge(u_gfix, g, n_gf, params.param.j_decay, params.param.GFAccu, params.param.GFMax,
		params.param.OrDo, params.param. OrPara);
    
      // Write out what is done
      push(xml_out,"Gauge_fixing_parameters");
      write(xml_out, "GFAccu",params.param.GFAccu);
      write(xml_out, "GFMax",params.param.GFMax);
      write(xml_out, "iterations",n_gf);
      pop(xml_out);
  
      // Check if the smeared gauge field is unitary
      unitarityCheck(u_gfix);
  
      // Again calculate some gauge invariant observables
      MesPlq(xml_out, "Gfix_observables", u_gfix);

      // Now store the configuration to a memory object
      {
	XMLBufferWriter file_xml, record_xml;
	push(file_xml, "gauge");
	write(file_xml, "id", int(0));
	pop(file_xml);
	record_xml << gauge_xml;

	// Store the gauge field
	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.gfix_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gfix_id) = u_gfix;
	TheNamedObjMap::Instance().get(params.named_obj.gfix_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.gfix_id).setRecordXML(record_xml);

	// Store the gauge rotation fields
	TheNamedObjMap::Instance().create< LatticeColorMatrix >(params.named_obj.gauge_rot_id);
	TheNamedObjMap::Instance().getData< LatticeColorMatrix >(params.named_obj.gauge_rot_id) = g;
	TheNamedObjMap::Instance().get(params.named_obj.gauge_rot_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_rot_id).setRecordXML(record_xml);
      }

      pop(xml_out);

      snoop.stop();
      QDPIO::cout << InlineCoulGaugeEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineCoulGaugeEnv::name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
