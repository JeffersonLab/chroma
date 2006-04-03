// $Id: inline_hyp_smear.cc,v 3.0 2006-04-03 04:59:04 edwards Exp $
/*! \file
 *  \brief Inline Hyp smearing
 */

#include "meas/inline/smear/inline_hyp_smear.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/hyp_smear.h"
#include "meas/smear/hyp_smear3d.h"
#include "util/info/proginfo.h"
#include "util/gauge/unit_check.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/default_gauge_field.h"

#include <sys/time.h>   // for timings

namespace Chroma 
{ 
  namespace InlineHypSmearEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineHypSmear(InlineHypSmearParams(xml_in, path));
    }

    const std::string name = "HYP_SMEAR";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Target file
  void read(XMLReader& xml, const string& path, InlineHypSmearParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    input.gauge_id = InlineDefaultGaugeField::readGaugeId(inputtop, "gauge_id");
    read(inputtop, "hyp_id", input.hyp_id);
  }


  //! Target file
  void write(XMLWriter& xml, const string& path, const InlineHypSmearParams::NamedObject_t& input)
  {
    push(xml, path);
    write(xml, "gauge_id", input.gauge_id);
    write(xml, "hyp_id", input.hyp_id);
    pop(xml);
  }


  //! Parameters for running code
  void read(XMLReader& xml, const string& path, InlineHypSmearParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);
    param.num_smear = 1;
    param.j_decay = -1;

    switch (version) 
    {
    case 2:
      break;

    case 3:
      read(paramtop, "num_smear", param.num_smear);
      break;

    case 4:
      read(paramtop, "num_smear", param.num_smear);
      read(paramtop, "j_decay", param.j_decay);
      break;

    default :
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "alpha1", param.alpha1);
    read(paramtop, "alpha2", param.alpha2);
    read(paramtop, "alpha3", param.alpha3);

    read(paramtop, "nrow", param.nrow);
  }

  //! Parameters
  void write(XMLWriter& xml, const string& path, const InlineHypSmearParams::Param_t& param)
  {
    push(xml, path);

    int version = 4;
    write(xml, "version", version);

    /* this version allows a variable num_smear */
    write(xml, "num_smear", param.num_smear);
    write(xml, "alpha1", param.alpha1);
    write(xml, "alpha2", param.alpha2);
    write(xml, "alpha3", param.alpha3);
    write(xml, "j_decay", param.num_smear);

    write(xml, "nrow", param.nrow);

    pop(xml);
  }



  // Param stuff
  InlineHypSmearParams::InlineHypSmearParams()
  { 
    frequency = 0; 
    named_obj.gauge_id = InlineDefaultGaugeField::getId();
  }

  InlineHypSmearParams::InlineHypSmearParams(XMLReader& xml_in, const std::string& path) 
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

      // Ids
      read(paramtop, "NamedObject", named_obj);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << InlineHypSmearEnv::name << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineHypSmearParams::write(XMLWriter& xml, const std::string& path) 
  {
    push(xml, path);

    Chroma::write(xml, "Param", param);
    Chroma::write(xml, "NamedObject", named_obj);

    pop(xml);
  }


  void 
  InlineHypSmear::operator()(unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Grab the gauge field
    XMLBufferWriter gauge_xml;
    multi1d<LatticeColorMatrix> u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
    TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

    push(xml_out, "hypsmear");
    write(xml_out, "update_no", update_no);
    
    QDPIO::cout << InlineHypSmearEnv::name << ": HYP smearing of gauge config" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);


    // Check if the gauge field configuration is unitarized
    clock_t t1 = clock();
    unitarityCheck(u);
    clock_t t2 = clock();
    QDPIO::cout << "Unitarity took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;
  

    // Calculate some gauge invariant observables just for info.
    t1 = clock();
    MesPlq(xml_out, "Observables", u);
    t2 = clock();
    QDPIO::cout << "Plaquette took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;


    // Now hyp smear
    multi1d<LatticeColorMatrix> u_hyp(Nd);

    Real BlkAccu = 1.0e-5;
    int BlkMax = 100;

    t1 = clock();
    if (params.param.num_smear > 0)
    {
      for (int n = 0; n < params.param.num_smear; n++)
      {
	if (params.param.j_decay < 0 || params.param.j_decay >= Nd)
	  Hyp_Smear(u, u_hyp, 
		    params.param.alpha1, params.param.alpha2, params.param.alpha3, 
		    BlkAccu, BlkMax);
	else
	  Hyp_Smear3d(u, u_hyp, 
		      params.param.alpha1, params.param.alpha2, params.param.alpha3, 
		      BlkAccu, BlkMax, params.param.j_decay);
      }
    }
    else
    {
      u_hyp = u;
    }
    t2 = clock();
    QDPIO::cout << "Hypsmear took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;

    // Calculate some gauge invariant observables just for info.

    MesPlq(xml_out, "HYP_observables", u_hyp);

    // Now store the configuration to a memory object
    {
      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      record_xml << gauge_xml;

      // Store the gauge field
      TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.hyp_id);
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.hyp_id) = u_hyp;
      TheNamedObjMap::Instance().get(params.named_obj.hyp_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.hyp_id).setRecordXML(record_xml);
    }
    pop(xml_out);

    snoop.stop();
    QDPIO::cout << InlineHypSmearEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineHypSmearEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
