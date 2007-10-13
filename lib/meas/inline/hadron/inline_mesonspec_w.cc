// $Id: inline_mesonspec_w.cc,v 3.16 2007-10-13 20:46:29 edwards Exp $
/*! \file
 * \brief Inline construction of meson spectrum
 *
 * Meson spectrum calculations
 */

#include "handle.h"
#include "meas/inline/hadron/inline_mesonspec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"

#include "meas/hadron/spin_insertion_factory.h"
#include "meas/hadron/spin_insertion_aggregate.h"

namespace Chroma 
{ 
  namespace InlineMesonSpecEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMesonSpec(InlineMesonSpecParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MESON_SPECTRUM";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= SpinInsertionEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  //! Reader for parameters
  void read(XMLReader& xml, const string& path, 
	    InlineMesonSpecParams::Param_t& param)
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
  void write(XMLWriter& xml, const string& path, 
	     const InlineMesonSpecParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "mom2_max", param.mom2_max);
    write(xml, "avg_equiv_mom", param.avg_equiv_mom);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, 
	    InlineMesonSpecParams::NamedObject_t::Correlators_t::CorrelatorTerms_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "first_id", input.first_id);
    read(inputtop, "second_id", input.second_id);
    read(inputtop, "factor", input.factor);

    input.source_spin_insertion = readXMLGroup(inputtop, "SourceSpinInsertion", "SpinInsertionType");
    input.sink_spin_insertion   = readXMLGroup(inputtop, "SinkSpinInsertion", "SpinInsertionType");
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, 
	     const InlineMesonSpecParams::NamedObject_t::Correlators_t::CorrelatorTerms_t& input)
  {
    push(xml, path);

    write(xml, "first_id", input.first_id);
    write(xml, "second_id", input.second_id);
    xml << input.source_spin_insertion.xml;
    xml << input.sink_spin_insertion.xml;
    write(xml, "factor", input.factor);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, 
	    InlineMesonSpecParams::NamedObject_t::Correlators_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "source_particle", input.source_particle);
    read(inputtop, "source_wavetype", input.source_wavetype);
    read(inputtop, "sink_particle", input.sink_particle);
    read(inputtop, "sink_wavetype", input.sink_wavetype);
    read(inputtop, "correlator_terms", input.correlator_terms);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, 
	     const InlineMesonSpecParams::NamedObject_t::Correlators_t& input)
  {
    push(xml, path);

    write(xml, "source_particle", input.source_particle);
    write(xml, "source_wavetype", input.source_wavetype);
    write(xml, "sink_particle", input.sink_particle);
    write(xml, "sink_wavetype", input.sink_wavetype);
    write(xml, "correlator_terms", input.correlator_terms);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, 
	    InlineMesonSpecParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "correlators", input.correlators);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, 
	     const InlineMesonSpecParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "correlators", input.correlators);

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


  // Writer
  void
  InlineMesonSpecParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Anonymous namespace
  namespace 
  {
    //! Useful structure holding sink props
    struct SinkPropContainer_t
    {
      ForwardProp_t prop_header;
      string quark_propagator_id;
      Real Mass;
    
      // Now loop over the various fermion masses
      string source_type;
      string source_disp_type;
      string sink_type;
      string sink_disp_type;
    };


    //! Useful structure holding sink props
    struct AllCorrelatorTerms_t
    {
      SinkPropContainer_t  sink_prop_1;
      SinkPropContainer_t  sink_prop_2;
    };


    //! Read a sink prop
    void readSinkProp(SinkPropContainer_t& s, const std::string& id)
    {
      try
      {
	// Try a cast to see if it succeeds
	const LatticePropagator& foo = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(id);

	// Snarf the data into a copy
	s.quark_propagator_id = id;
	
	// Snarf the prop info. This is will throw if the prop_id is not there
	XMLReader prop_file_xml, prop_record_xml;
	TheNamedObjMap::Instance().get(id).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(id).getRecordXML(prop_record_xml);
   
	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	{
	  read(prop_record_xml, "/SinkSmear", s.prop_header);
	  
	  read(prop_record_xml, "/SinkSmear/PropSource/Source/SourceType", s.source_type);
	  read(prop_record_xml, "/SinkSmear/PropSource/Source/Displacement/DisplacementType", 
	       s.source_disp_type);

	  read(prop_record_xml, "/SinkSmear/PropSink/Sink/SinkType", s.sink_type);
	  read(prop_record_xml, "/SinkSmear/PropSink/Sink/Displacement/DisplacementType", 
	       s.sink_disp_type);
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


      // Derived from input prop
      // Hunt around to find the mass
      // NOTE: this may be problematic in the future if actions are used with no
      // clear def. of a Mass
      QDPIO::cout << "Try action and mass" << endl;
      s.Mass = getMass(s.prop_header.prop_header.fermact);
    
      QDPIO::cout << "FermAct = " << s.prop_header.prop_header.fermact.id << endl;
      QDPIO::cout << "Mass = " << s.Mass << endl;
    }


    //! Read all sinks
    void readAllSinks(multi1d<AllCorrelatorTerms_t>& s, 
		      const multi1d<InlineMesonSpecParams::NamedObject_t::Correlators_t::CorrelatorTerms_t>& correlator_terms)
    {
      s.resize(correlator_terms.size());

      for(int i=0; i < correlator_terms.size(); ++i)
      {
	QDPIO::cout << "Attempt to parse forward propagator = " << correlator_terms[i].first_id << endl;
	readSinkProp(s[i].sink_prop_1, correlator_terms[i].first_id);
	QDPIO::cout << "Forward propagator successfully parsed" << endl;

	QDPIO::cout << "Attempt to parse forward propagator = " << correlator_terms[i].second_id << endl;
	readSinkProp(s[i].sink_prop_2, correlator_terms[i].second_id);
	QDPIO::cout << "Forward propagator successfully parsed" << endl;
      }
    }

  } // namespace anonymous



  // Function call
  void 
  InlineMesonSpec::operator()(unsigned long update_no,
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
      func(update_no, xml);
    }
    else
    {
      func(update_no, xml_out);
    }
  }



  // Real work done here
  void 
  InlineMesonSpec::func(unsigned long update_no,
			XMLWriter& xml_out) 
  {
    START_CODE();

    QDPIO::cout << InlineMesonSpecEnv::name << ": meson spectroscopy for Wilson-like fermions" << endl;

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineMesonSpecEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMesonSpecEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "MesonSpectrum");
    write(xml_out, "update_no", update_no);

    proginfo(xml_out);    // Print out basic program info

    QDPIO::cout << "xml_file = " << params.xml_file << endl;

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 3);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    push(xml_out, "Wilson_hadron_measurements");

    // Now loop over the various fermion pairs
    for(int lpair=0; lpair < params.named_obj.correlators.size(); ++lpair)
    {
      const InlineMesonSpecParams::NamedObject_t::Correlators_t named_obj = params.named_obj.correlators[lpair];

      push(xml_out, "elem");

      write(xml_out, "source_particle", named_obj.source_particle);
      write(xml_out, "source_wavetype", named_obj.source_wavetype);
      write(xml_out, "sink_particle", named_obj.sink_particle);
      write(xml_out, "sink_wavetype", named_obj.sink_wavetype);

      multi1d<AllCorrelatorTerms_t> all_sinks;
      readAllSinks(all_sinks, named_obj.correlator_terms);

      // Derived from input prop
      multi1d<int> t_srce
                  = all_sinks[0].sink_prop_1.prop_header.source_header.getTSrce();
      int j_decay = all_sinks[0].sink_prop_1.prop_header.source_header.j_decay;
      int t0      = all_sinks[0].sink_prop_1.prop_header.source_header.t_source;

      // Sanity checks
      for(int loop=0; loop < all_sinks.size(); ++loop)
      {
	if (all_sinks[loop].sink_prop_2.prop_header.source_header.j_decay != j_decay)
	{
	  QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
	if (all_sinks[loop].sink_prop_2.prop_header.source_header.t_source != 
	    all_sinks[loop].sink_prop_1.prop_header.source_header.t_source)
	{
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
      }
	

      // Initialize the slow Fourier transform phases
      SftMom phases(params.param.mom2_max, t_srce, params.param.avg_equiv_mom,
                    j_decay);

      // Keep a copy of the phases with NO momenta
      SftMom phases_nomom(0, true, j_decay);

      // Masses
      write(xml_out, "Mass_1", all_sinks[0].sink_prop_1.Mass);
      write(xml_out, "Mass_2", all_sinks[0].sink_prop_2.Mass);
      write(xml_out, "t0", t0);

      // Save prop input
      push(xml_out, "Forward_prop_headers");
      for(int loop=0; loop < all_sinks.size(); ++loop)
      {
	push(xml_out, "elem");
	write(xml_out, "First_forward_prop", all_sinks[loop].sink_prop_1.prop_header);
	write(xml_out, "Second_forward_prop", all_sinks[loop].sink_prop_2.prop_header);
	pop(xml_out);
      }
      pop(xml_out);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      push(xml_out, "Forward_prop_correlator");
      for(int loop=0; loop < all_sinks.size(); ++loop)
      {
	const LatticePropagator& sink_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks[loop].sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks[loop].sink_prop_2.quark_propagator_id);

	push(xml_out, "elem");
	write(xml_out, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	write(xml_out, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
	pop(xml_out);
      }
      pop(xml_out);


      push(xml_out, "SourceSinkType");
      for(int loop=0; loop < all_sinks.size(); ++loop)
      {
	push(xml_out, "elem");
	QDPIO::cout << "Source_type_1 = " << all_sinks[loop].sink_prop_1.source_type << endl;
	QDPIO::cout << "Sink_type_1 = " << all_sinks[loop].sink_prop_1.sink_type << endl;
	QDPIO::cout << "Source_type_2 = " << all_sinks[loop].sink_prop_2.source_type << endl;
	QDPIO::cout << "Sink_type_2 = " << all_sinks[loop].sink_prop_2.sink_type << endl;

	write(xml_out, "source_type_1", all_sinks[loop].sink_prop_1.source_type);
	write(xml_out, "source_disp_type_1", all_sinks[loop].sink_prop_1.source_disp_type);
	write(xml_out, "sink_type_1", all_sinks[loop].sink_prop_1.sink_type);
	write(xml_out, "sink_disp_type_1", all_sinks[loop].sink_prop_1.sink_disp_type);

	write(xml_out, "source_type_2", all_sinks[loop].sink_prop_2.source_type);
	write(xml_out, "source_disp_type_2", all_sinks[loop].sink_prop_2.source_disp_type);
	write(xml_out, "sink_type_2", all_sinks[loop].sink_prop_2.sink_type);
	write(xml_out, "sink_disp_type_2", all_sinks[loop].sink_prop_2.sink_disp_type);
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
	LatticeComplex corr_fn = zero;

	for(int loop=0; loop < all_sinks.size(); ++loop)
	{
	  const LatticePropagator& sink_prop_1 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks[loop].sink_prop_1.quark_propagator_id);
	  const LatticePropagator& sink_prop_2 = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks[loop].sink_prop_2.quark_propagator_id);

	  LatticePropagator prop_1, prop_2;

	  // Factory constructions
	  try
	  {
	    // Create the source spin insertion object
	    {
	      std::istringstream  xml_s(named_obj.correlator_terms[loop].source_spin_insertion.xml);
	      XMLReader  inserttop(xml_s);
//	      const string insert_path = "/SourceSpinInsertion";
	
	      Handle< SpinInsertion<LatticePropagator> > sourceSpinInsertion(
		ThePropSpinInsertionFactory::Instance().createObject(
		  named_obj.correlator_terms[loop].source_spin_insertion.id,
		  inserttop,
		  named_obj.correlator_terms[loop].source_spin_insertion.path));

	      prop_1 = (*sourceSpinInsertion)(sink_prop_1);
	    }

	    // Create the sink spin insertion object
	    {
	      std::istringstream  xml_s(named_obj.correlator_terms[loop].sink_spin_insertion.xml);
	      XMLReader  inserttop(xml_s);
//	      const string insert_path = "/SinkSpinInsertion";
	
	      Handle< SpinInsertion<LatticePropagator> > sinkSpinInsertion(
		ThePropSpinInsertionFactory::Instance().createObject(
		  named_obj.correlator_terms[loop].sink_spin_insertion.id,
		  inserttop,
		  named_obj.correlator_terms[loop].sink_spin_insertion.path));

	      LatticePropagator prop_tmp = Gamma(G5) * adj(sink_prop_2) * Gamma(G5);

	      prop_2 = (*sinkSpinInsertion)(prop_tmp);
	    }
	  }
	  catch(const std::string& e) 
	  {
	    QDPIO::cerr << InlineMesonSpecEnv::name << ": Caught Exception inserting: " 
			<< e << endl;
	    QDP_abort(1);
	  }


	  LatticeComplex tmp = trace(prop_2 * prop_1);

	  corr_fn += named_obj.correlator_terms[loop].factor * tmp;
	}

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

	  multi1d<DComplex> mesprop(length);
	  for (int t=0; t < length; ++t) 
	  {
	    int t_eff = (t - t0 + length) % length;
	    mesprop[t_eff] = hsum[sink_mom_num][t];
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
