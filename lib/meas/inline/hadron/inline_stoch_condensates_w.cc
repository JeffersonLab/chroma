// $Id: inline_stoch_condensates_w.cc,v 3.1 2007-08-10 21:27:02 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic condensates
 *
 */

#include "handle.h"
#include "meas/inline/hadron/inline_stoch_condensates_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineStochCondensatesEnv 
  { 
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

    const std::string name = "STOCH_CONDENSATES";

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


    // Reader for input parameters
    void read(XMLReader& xml, const string& path, InlineStochCondensatesEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	/**************************************************************************/
	break;

      default :
	/**************************************************************************/

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", param.mom2_max);
    }


    // Reader for input parameters
    void write(XMLWriter& xml, const string& path, const InlineStochCondensatesEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);

      pop(xml);
    }


    //! Propagator parameters
    void read(XMLReader& xml, const string& path, InlineStochCondensatesEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "soln_files", input.soln_files);
    }

    //! Propagator parameters
    void write(XMLWriter& xml, const string& path, const InlineStochCondensatesEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "soln_files", input.soln_files);

      pop(xml);
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
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "Param", param);

      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }



    //--------------------------------------------------------------

    //! Structure holding a source and its solutions
    struct QuarkSourceSolutions_t
    {
      //! Structure holding solutions
      struct QuarkSolution_t
      {
	LatticeFermion     source;
	LatticeFermion     soln;
	PropSourceConst_t  source_header;
	ChromaProp_t       prop_header;
      };

      int   j_decay;
      Seed  seed;
      multi1d<QuarkSolution_t>  dilutions;
    };


    //--------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "stoch_meson");
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


    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      StopWatch swatch;

      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineStochCondensatesEnv::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << InlineStochCondensatesEnv::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "stoch_condensates");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineStochCondensatesEnv::name << ": Stochastic Condensates" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // First calculate some gauge invariant observables just for info.
      // This is really cheap.
      MesPlq(xml_out, "Observables", u);

      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);

      //
      // Read the source and solutions
      //
      swatch.reset();
      swatch.start();
      QuarkSourceSolutions_t  quark;

      try
      {
	QDPIO::cout << "Attempt to read solutions" << endl;
	quark.dilutions.resize(params.named_obj.soln_files.size());

	QDPIO::cout << "dilutions.size= " << quark.dilutions.size() << endl;
	for(int i=0; i < quark.dilutions.size(); ++i)
	{
	  XMLReader file_xml, record_xml;

	  QDPIO::cout << "reading file= " << params.named_obj.soln_files[i] << endl;
	  QDPFileReader from(file_xml, params.named_obj.soln_files[i], QDPIO_SERIAL);
	  read(from, record_xml, quark.dilutions[i].soln);
	  close(from);
	
	  read(record_xml, "/Propagator/PropSource", quark.dilutions[i].source_header);
	  read(record_xml, "/Propagator/ForwardProp", quark.dilutions[i].prop_header);
	}
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting headers: " << e << endl;
	QDP_abort(1);
      }
      swatch.stop();

      QDPIO::cout << "Sources and solutions successfully read: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;



      //
      // Check for each quark source that the solutions have their diluted
      // on every site only once
      //
      swatch.reset();
      swatch.start();

      try
      {
	push(xml_out, "Norms");
	bool first = true;
	int  N;
	LatticeFermion quark_noise;      // noisy source on entire lattice

	for(int i=0; i < quark.dilutions.size(); ++i)
	{
	  std::istringstream  xml_s(quark.dilutions[i].source_header.source.xml);
	  XMLReader  sourcetop(xml_s);
//	QDPIO::cout << "Source = " << quark.dilutions[i].source_header.source.id << endl;

	  if (quark.dilutions[i].source_header.source.id != DiluteZNQuarkSourceConstEnv::name)
	  {
	    QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
	    QDP_abort(1);
	  }

	  QDPIO::cout << "Dilution num= " << i << endl;

	  // Manually create the params so I can peek into them and use the source constructor
	  DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
							 quark.dilutions[i].source_header.source.path);
	  DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
	  if (first) 
	  {
	    first = false;

	    quark.j_decay = srcParams.j_decay;

	    // Grab N
	    N = srcParams.N;

	    // Set the seed to desired value
	    quark.seed = srcParams.ran_seed;
	    QDP::RNG::setrn(quark.seed);

	    // Create the noisy quark source on the entire lattice
	    zN_src(quark_noise, N);
	  }

	  // The seeds must always agree - here the seed is the unique id of the source
	  if ( toBool(srcParams.ran_seed != quark.seed) )
	  {
	    QDPIO::cerr << "dilution=" << i << " seed does not match" << endl;
	    QDP_abort(1);
	  }

	  // The N's must always agree
	  if ( toBool(srcParams.N != N) )
	  {
	    QDPIO::cerr << "dilution=" << i << " N does not match" << endl;
	    QDP_abort(1);
	  }

	  // Use a trick here, create the source and subtract it from the global noisy
	  // Check at the end that the global noisy is zero everywhere.
	  // NOTE: the seed will be set every call
	  quark.dilutions[i].source = srcConst(u);
	  quark_noise -= quark.dilutions[i].source;

#if 0
	  // Diagnostic
	  {
	    // Keep a copy of the phases with NO momenta
	    SftMom phases_nomom(0, true, quark.dilutions[i].source_header.j_decay);

	      multi1d<Double> source_corr = sumMulti(localNorm2(quark.dilutions[i].source), 
						     phases_nomom.getSet());
	      
	      multi1d<Double> soln_corr = sumMulti(localNorm2(quark.dilutions[i].soln), 
						   phases_nomom.getSet());

	      push(xml_out, "elem");
	      write(xml_out, "n", n);
	      write(xml_out, "i", i);
	      write(xml_out, "source_corr", source_corr);
	      write(xml_out, "soln_corr", soln_corr);
	      pop(xml_out);
	    }
#endif
	  } // end for i

	  Double dcnt = norm2(quark_noise);
	  if (toDouble(dcnt) != 0.0)  // problematic - seems to work with unnormalized sources 
	  {
	    QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	    QDP_abort(1);
	  }

	pop(xml_out);
      } // end try
      catch(const std::string& e) 
      {
	QDPIO::cerr << ": Caught Exception creating source: " << e << endl;
	QDP_abort(1);
      }

      swatch.stop();

      QDPIO::cout << "Sources saturated: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


      //
      // Meson operators
      //
      int j_decay = quark.j_decay;

      // Initialize the slow Fourier transform phases
      SftMom phases(params.param.mom2_max, false, j_decay);
    
      // Length of lattice in decay direction
      int length = phases.numSubsets();

      // Start operator contractions
      swatch.reset();
      swatch.start();

      push(xml_out, "ContractionParams");
      write(xml_out, "mom2_max", params.param.mom2_max);
      write(xml_out, "decay_dir", j_decay);
      pop(xml_out);

      // Construct operator A
      int G5 = Ns*Ns-1;

      push(xml_out, "SimpleGamma");

      for(int gamma_value=0; gamma_value < Ns*Ns; ++gamma_value)
      {
	push(xml_out, "elem");
	write(xml_out, "gamma_value", gamma_value);

	LatticeComplex corr_fn = zero;

	for(int i=0; i < quark.dilutions.size(); ++i)
	{
	  corr_fn += 
	    localInnerProduct(quark.dilutions[i].source, Gamma(gamma_value) * quark.dilutions[i].soln);
	} // end for i

	multi2d<DComplex> hsum;
	hsum = phases.sft(corr_fn);

	// Loop over sink momenta
	push(xml_out, "momenta");

	for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
	{
	  push(xml_out, "elem");
	  write(xml_out, "sink_mom_num", sink_mom_num);
	  write(xml_out, "sink_mom", phases.numToMom(sink_mom_num));

	  multi1d<DComplex> corr(length);
	  for (int t=0; t < length; ++t) 
	  {
	    corr[t] = hsum[sink_mom_num][t];
	  }

	  write(xml_out, "condensate", corr);
	  pop(xml_out);
	  
	} // end for(sink_mom_num)
 
	pop(xml_out); // momenta
	pop(xml_out); // elem
      } // end for g
      
      pop(xml_out); // SimpleGamma

      swatch.stop();

      QDPIO::cout << "Simple gamma condensates computed: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;

      // Close the namelist output file XMLDAT
      pop(xml_out);     // StochCondensates

      snoop.stop();
      QDPIO::cout << InlineStochCondensatesEnv::name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << InlineStochCondensatesEnv::name << ": ran successfully" << endl;

      END_CODE();
    } 

  } // namespace InlineStochCondensatesEnv

} // namespace Chroma
