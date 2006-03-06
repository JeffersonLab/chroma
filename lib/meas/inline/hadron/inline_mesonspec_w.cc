// $Id: inline_mesonspec_w.cc,v 2.6 2006-03-06 22:14:11 edwards Exp $
/*! \file
 * \brief Inline construction of meson spectrum
 *
 * Meson spectrum calculations
 */

#include "meas/inline/hadron/inline_mesonspec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineMesonSpecEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineMesonSpec(InlineMesonSpecParams(xml_in, path));
    }

    const std::string name = "MESON_SPECTRUM";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };



  //! Reader for parameters
  void read(XMLReader& xml, const string& path, InlineMesonSpecParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      break;

    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "mom2_max", param.mom2_max);
    read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);
  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineMesonSpecParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "mom2_max", param.mom2_max);
    write(xml, "avg_equiv_mom", param.avg_equiv_mom);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineMesonSpecParams::NamedObject_t::Props_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "sink_ids", input.sink_ids);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineMesonSpecParams::NamedObject_t::Props_t& input)
  {
    push(xml, path);

    write(xml, "sink_ids", input.sink_ids);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineMesonSpecParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_ids", input.prop_ids);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineMesonSpecParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "prop_ids", input.prop_ids);

    pop(xml);
  }


  // Param stuff
  InlineMesonSpecParams::InlineMesonSpecParams() { frequency = 0; }

  InlineMesonSpecParams::InlineMesonSpecParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineMesonSpecParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Function call
  void 
  InlineMesonSpec::operator()(const multi1d<LatticeColorMatrix>& u,
			      XMLBufferWriter& gauge_xml,
			      unsigned long update_no,
			      XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "MesonSpectrum");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(u, gauge_xml, update_no, xml);
    }
    else
    {
      func(u, gauge_xml, update_no, xml_out);
    }
  }


  // Real work done here
  void 
  InlineMesonSpec::func(const multi1d<LatticeColorMatrix>& u,
			XMLBufferWriter& gauge_xml,
			unsigned long update_no,
			XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    push(xml_out, "MesonSpectrum");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " MESONSPEC: Meson spectroscopy for Wilson-like fermions" << endl;
    QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;
    QDPIO::cout << "     volume: " << QDP::Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << QDP::Layout::lattSize()[i];
    }
    QDPIO::cout << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    push(xml_out, "Wilson_hadron_measurements");

    // Now loop over the various fermion masses
    for (int lpair=0; lpair < params.named_obj.prop_ids.size(); ++lpair)
    {
      const InlineMesonSpecParams::NamedObject_t::Props_t named_obj = params.named_obj.prop_ids[lpair];

      push(xml_out, "elem");

      // Require 2 props
      const int Nprops = 2;
      if (named_obj.sink_ids.size() != Nprops)
      {
	QDPIO::cerr << InlineMesonSpecEnv::name << ": needs a 2 element array of sink propagators" << endl;
	QDP_abort(2);
      }

      multi1d<ForwardProp_t> prop_header(Nprops);
      multi1d<LatticePropagator> quark_propagator(Nprops);
      multi1d<Real> Mass(Nprops);
    
      multi2d<int> bc(Nprops, 4); 
    
      // Now loop over the various fermion masses
      multi1d<string> source_type(Nprops);
      multi1d<string> source_disp_type(Nprops);
      multi1d<string> sink_type(Nprops);
      multi1d<string> sink_disp_type(Nprops);

      for(int loop=0; loop < Nprops; ++loop)
      {
	QDPIO::cout << "Attempt to parse forward propagator = " << named_obj.sink_ids[loop] << endl;
	try
	{
	  // Snarf the data into a copy
	  quark_propagator[loop] =
	    TheNamedObjMap::Instance().getData<LatticePropagator>(named_obj.sink_ids[loop]);
	
	  // Snarf the prop info. This is will throw if the prop_id is not there
	  XMLReader prop_file_xml, prop_record_xml;
	  TheNamedObjMap::Instance().get(named_obj.sink_ids[loop]).getFileXML(prop_file_xml);
	  TheNamedObjMap::Instance().get(named_obj.sink_ids[loop]).getRecordXML(prop_record_xml);
   
	  // Try to invert this record XML into a ChromaProp struct
	  // Also pull out the id of this source
	  {
	    read(prop_record_xml, "/SinkSmear", prop_header[loop]);

	    read(prop_record_xml, "/SinkSmear/PropSource/Source/SourceType", source_type[loop]);
	    read(prop_record_xml, "/SinkSmear/PropSource/Source/Displacement/DisplacementType", 
		 source_disp_type[loop]);

	    read(prop_record_xml, "/SinkSmear/PropSink/Sink/SinkType", sink_type[loop]);
	    read(prop_record_xml, "/SinkSmear/PropSink/Sink/Displacement/DisplacementType", 
		 sink_disp_type[loop]);
	  }
	}
	catch( std::bad_cast ) 
	{
	  QDPIO::cerr << InlineMesonSpecEnv::name << ": caught dynamic cast error" 
		      << endl;
	  QDP_abort(1);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << InlineMesonSpecEnv::name << ": error message: " << e 
		      << endl;
	  QDP_abort(1);
	}
	QDPIO::cout << "Forward propagator successfully parsed" << endl;


	// Derived from input prop
	// Hunt around to find the mass
	// NOTE: this may be problematic in the future if actions are used with no
	// clear def. of a Mass
	std::istringstream  xml_s(prop_header[loop].prop_header.fermact);
	XMLReader  fermacttop(xml_s);
	const string fermact_path = "/FermionAction";
	string fermact;
      
	QDPIO::cout << "Try action and mass" << endl;
	try
	{
	  XMLReader top(fermacttop, fermact_path);

	  read(top, "FermAct", fermact);

	  // Yuk - need to hop some hoops. This should be isolated.
	  if (top.count("Mass") != 0) 
	  {
	    read(top, "Mass", Mass[loop]);
	  }
	  else if (top.count("Kappa") != 0)
	  {
	    Real Kappa;
	    read(top, "Kappa", Kappa);
	    Mass[loop] = kappaToMass(Kappa);    // Convert Kappa to Mass
	  }
	  else if (top.count("m_q") != 0) 
	  {
	    read(top, "m_q", Mass[loop]);
	  }
	  else
	  {
	    QDPIO::cerr << "Neither Mass nor Kappa found" << endl;
	    throw std::string("Neither Mass nor Kappa found");
	  }
	  bc[loop] = getFermActBoundary(prop_header[loop].prop_header.fermact);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << "Error reading fermact or mass: " << e << endl;
	  QDP_abort(1);
	}
    
	QDPIO::cout << "FermAct = " << fermact << endl;
	QDPIO::cout << "Mass = " << Mass[loop] << endl;
      }

 
      //
      // Sanity checks
      // For now, force the type of source and sink smearings to agree.
      // This could/should be relaxed to say point/smeared for BNDST
      // type constructions in which case the BNDST type will be used
      // in place of point or shell.
      //
      // NOTE: the only real requirement is that in BARHQLQ that the
      // width of sink smearing is the same since no antisymmetrization
      // is done.
      //
      for (int loop(1); loop < Nprops; ++loop)
      {
	if (source_type[loop] != source_type[loop])
	{
	  QDPIO::cerr << InlineMesonSpecEnv::name << ": prop source smearing types do not agree" << endl;
	  QDP_abort(1);
	}

	if (sink_type[loop] != sink_type[0])
	{
	  QDPIO::cerr << InlineMesonSpecEnv::name << ": prop sink smearing types do not agree" << endl;
	  QDP_abort(1);
	}
      }


      // Derived from input prop
      int j_decay = prop_header[0].source_header.j_decay;
      int t0      = prop_header[0].source_header.t_source;
      int bc_spec = bc[0][j_decay] ;
      for (int loop(0); loop < Nprops; ++loop)
      {
	if (prop_header[loop].source_header.j_decay != j_decay)
	{
	  QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
	if (bc[loop][j_decay] != bc_spec)
	{
	  QDPIO::cerr << "Error!! bc must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
	if (prop_header[loop].source_header.t_source != prop_header[0].source_header.t_source)
	{
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
      }
  

      // Initialize the slow Fourier transform phases
      SftMom phases(params.param.mom2_max, params.param.avg_equiv_mom, j_decay);

      // Keep a copy of the phases with NO momenta
      SftMom phases_nomom(0, true, j_decay);

      // Next array element - name auto-written
      write(xml_out, "Masses", Mass);
      write(xml_out, "t0", t0);

      // Save prop input
      write(xml_out, "ForwardProp", prop_header);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	multi1d< multi1d<Double> > forward_prop_corr(Nprops);
	for (int loop=0; loop < Nprops; ++loop)
	{
	  forward_prop_corr[loop] = sumMulti(localNorm2(quark_propagator[loop]), 
					     phases.getSet());
	}
	  
	push(xml_out, "Forward_prop_correlator");
	write(xml_out, "forward_prop_corr", forward_prop_corr);
	pop(xml_out);
      }


      // Construct group name for output
      string src_type;
      if (source_type[0] == "POINT_SOURCE")
	src_type = "Point";
      else if (source_type[0] == "SHELL_SOURCE")
	src_type = "Shell";
      else if (source_type[0] == "WALL_SOURCE")
	src_type = "Wall";
      else
      {
	QDPIO::cerr << "Unsupported source type" << endl;
	QDP_abort(1);
      }

      string snk_type;
      if (sink_type[0] == "POINT_SINK")
	snk_type = "Point";
      else if (sink_type[0] == "SHELL_SINK")
	snk_type = "Shell";
      else if (sink_type[0] == "WALL_SINK")
	snk_type = "Wall";
      else
      {
	QDPIO::cerr << "Unsupported sink type" << endl;
	QDP_abort(1);
      }

      string source_sink_type = src_type + "_" + snk_type;
      QDPIO::cout << "Source type = " << src_type << endl;
      QDPIO::cout << "Sink type = "   << snk_type << endl;

      push(xml_out, "SourceSinkType");
      for (int loop=0; loop < Nprops; ++loop)
      {
	push(xml_out, "elem");
	write(xml_out, "quark_number", loop);
	write(xml_out, "source_type", source_type[loop]);
	write(xml_out, "source_disp_type", source_disp_type[loop]);
	write(xml_out, "sink_type", sink_type[loop]);
	write(xml_out, "sink_disp_type", sink_disp_type[loop]);
	pop(xml_out);
      }
      pop(xml_out);

      // Do the mesons
      push(xml_out, "Mesons");
      {
	// Length of lattice in decay direction
	int length = phases.numSubsets();

	// Construct the anti-quark propagator from quark_propagator[1]
	int G5 = Ns*Ns-1;

	// Construct the meson correlation function
	LatticeComplex corr_fn =
	  trace((Gamma(G5) * adj(quark_propagator[1])  * Gamma(G5)) * quark_propagator[0]);

	multi2d<DComplex> hsum;
	hsum = phases.sft(corr_fn);

	// Loop over sink momenta
	XMLArrayWriter xml_sink_mom(xml_out,phases.numMom());
	push(xml_sink_mom, "momenta");

	for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
	{
	  push(xml_sink_mom);
	  write(xml_sink_mom, "sink_mom_num", sink_mom_num);
	  write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num));

	  multi1d<Real> mesprop(length);
	  for (int t=0; t < length; ++t) 
	  {
	    int t_eff = (t - t0 + length) % length;
	    mesprop[t_eff] = real(hsum[sink_mom_num][t]);
	  }

	  write(xml_sink_mom, "mesprop", mesprop);
	  pop(xml_sink_mom);
	
	} // end for(sink_mom_num)
 
	pop(xml_sink_mom);
      }
      pop(xml_out);  // Mesons

      pop(xml_out);  // array element
    }
    pop(xml_out);  // Wilson_spectroscopy
    pop(xml_out);  // mesonspec

    snoop.stop();
    QDPIO::cout << InlineMesonSpecEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineMesonSpecEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
