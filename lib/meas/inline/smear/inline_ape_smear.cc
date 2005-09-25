// $Id: inline_ape_smear.cc,v 2.0 2005-09-25 21:04:39 edwards Exp $
/*! \file
 *  \brief Inline APE smearing
 */

#include "meas/inline/smear/inline_ape_smear.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "util/info/proginfo.h"
#include "util/gauge/unit_check.h"
#include "meas/inline/io/named_objmap.h"


namespace Chroma 
{ 
  namespace InlineAPESmearEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineAPESmear(InlineAPESmearParams(xml_in, path));
    }

    const std::string name = "APE_SMEAR";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Parameters for running code
  void read(XMLReader& xml, const string& path, InlineAPESmearParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 2:
      break;

    default :

      QDPIO::cerr << "Input version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "link_smear_num", param.link_smear_num);
    read(paramtop, "link_smear_fact", param.link_smear_fact);
    read(paramtop, "j_decay", param.j_decay);
    read(paramtop, "nrow", param.nrow);
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const InlineAPESmearParams::Param_t& param)
  {
    push(xml, path);
    
    int version = 2;
    write(xml, "version", version);
    write(xml, "link_smear_num", param.link_smear_num);
    write(xml, "link_smear_fact", param.link_smear_fact);
    write(xml, "j_decay", param.j_decay);
    write(xml, "nrow", param.nrow);

    pop(xml);
  }

  // Reader for out gauge file
  void read(XMLReader& xml, const string& path, InlineAPESmearParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);
    read(inputtop, "ape_id", input.ape_id);
  }

  // Reader for out gauge file
  void write(XMLWriter& xml, const string& path, const InlineAPESmearParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "ape_id", input.ape_id);

    pop(xml);
  }


  // Param stuff
  InlineAPESmearParams::InlineAPESmearParams() { frequency = 0; }

  InlineAPESmearParams::InlineAPESmearParams(XMLReader& xml_in, const std::string& path) 
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

      // Read in the ape outfile
      read(paramtop, "NamedObject", named_obj);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  // Write params
  void
  InlineAPESmearParams::write(XMLWriter& xml, const std::string& path) 
  {
    push(xml, path);
      
    Chroma::write(xml, "Param", param);
    Chroma::write(xml, "NamedObject", named_obj);

    pop(xml);
  }


  void 
  InlineAPESmear::operator()(const multi1d<LatticeColorMatrix>& u,
			     XMLBufferWriter& gauge_xml,
			     unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    push(xml_out, "apesmear");
    write(xml_out, "update_no", update_no);
    
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


    // Now ape smear
    multi1d<LatticeColorMatrix> u_ape(Nd);
    u_ape = u;

    if (params.param.link_smear_num > 0)
    {
      QDPIO::cout << "APE Smear gauge field" << endl;

      int BlkMax = 100;
      Real BlkAccu = 1.0e-5;

      for(int i=0; i < params.param.link_smear_num; ++i)
      {
	multi1d<LatticeColorMatrix> u_tmp(Nd);

	for(int mu = 0; mu < Nd; ++mu)
	  if ( mu != params.param.j_decay )
	    APE_Smear(u_ape, u_tmp[mu], mu, 0,
		      params.param.link_smear_fact, BlkAccu, BlkMax,
		      params.param.j_decay);
	  else
	    u_tmp[mu] = u_ape[mu];

	u_ape = u_tmp;
      }
      QDPIO::cout << "Gauge field APE-smeared!" << endl;
    }

    // Write out what is done
    push(xml_out,"Smearing_parameters");
    write(xml_out, "link_smear_num",params.param.link_smear_num);
    write(xml_out, "link_smear_fact",params.param.link_smear_fact);
    write(xml_out, "j_decay",params.param.j_decay);
    pop(xml_out);
  
    // Check if the smeared gauge field is unitary
    unitarityCheck(u_ape);
  
    // Again calculate some gauge invariant observables
    MesPlq(xml_out, "Observables", u_ape);

    // Now store the configuration to a memory object
    {
      XMLBufferWriter file_xml, record_xml;
      push(file_xml, "gauge");
      write(file_xml, "id", int(0));
      pop(file_xml);
      record_xml << gauge_xml;

      // Store the gauge field
      TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.ape_id);
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.ape_id) = u_ape;
      TheNamedObjMap::Instance().get(params.named_obj.ape_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.ape_id).setRecordXML(record_xml);
    }

    pop(xml_out);

    END_CODE();
  } 

};
