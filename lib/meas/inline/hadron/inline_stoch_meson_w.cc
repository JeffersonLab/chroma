// $Id: inline_stoch_meson_w.cc,v 3.5 2009-08-05 15:26:49 jbulava Exp $
/*! \file
 * \brief Inline measurement of stochastic meson operator
 *
 */

#include "handle.h"
#include "meas/inline/hadron/inline_stoch_meson_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/sources/source_smearing_aggregate.h"
#include "meas/sources/source_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/smear/quark_source_sink.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineStochMesonEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineStochMeson(InlineStochMesonParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "STOCH_MESON";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= QuarkSourceSmearingEnv::registerAll();
	success &= QuarkSinkSmearingEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  // Operator parameters
  void read(XMLReader& xml, const string& path, InlineStochMesonParams::Prop_t::Operator_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "soln_files", input.soln_files);
  }


  // Operator parameters
  void write(XMLWriter& xml, const string& path, const InlineStochMesonParams::Prop_t::Operator_t& input)
  {
    push(xml, path);
    write(xml, "soln_files", input.soln_files);
    pop(xml);
  }


  // Propagator parameters
  void read(XMLReader& xml, const string& path, InlineStochMesonParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "operator_file", input.op_file);
    read(inputtop, "operator", input.op);
  }


  // Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineStochMesonParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "operator_file", input.op_file);
    write(xml, "operator", input.op);

    pop(xml);
  }


  // Reader for input parameters
  void read(XMLReader& xml, const string& path, InlineStochMesonParams::Param_t& param)
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
  void write(XMLWriter& xml, const string& path, const InlineStochMesonParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;

    write(xml, "version", version);
    write(xml, "mom2_max", param.mom2_max);

    pop(xml);
  }


  //! Propagator parameters
  void read(XMLReader& xml, const string& path, InlineStochMesonParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "Prop", input.prop);
  }

  //! Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineStochMesonParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "Prop", input.prop);

    pop(xml);
  }


  // Param stuff
  InlineStochMesonParams::InlineStochMesonParams()
  { 
    frequency = 0; 
  }

  InlineStochMesonParams::InlineStochMesonParams(XMLReader& xml_in, const std::string& path) 
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

      // Source smearing
      read(paramtop, "SourceSmearing", source_smearing);

      // Sink smearing
      read(paramtop, "SinkSmearing", sink_smearing);

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
  InlineStochMesonParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Source smearing
    Chroma::write(xml_out, "SourceSmearing", source_smearing);

    // Sink smearing
    Chroma::write(xml_out, "SinkSmearing", sink_smearing);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "NamedObject", named_obj);

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


  //! Meson operator
  struct MesonOperator_t
  {
    //! Serialize generalized operator object
    multi1d<Complex> serialize();

    //! Meson operator
    struct MesonOperatorInsertion_t
    {
      //! Meson operator element
      struct MesonOperatorElement_t
      {
	multi2d<DComplex> elem;              /*!< time slice and momenta number */
      };
    
      multi2d<MesonOperatorElement_t> op;    /*!< hybrid list indices */
    };

    std::string   smearing_l;        /*!< string holding smearing xml */
    std::string   smearing_r;        /*!< string holding smearing xml */

    int           mom2_max;
    int           j_decay;
    Seed          seed_l;
    Seed          seed_r;
    multi1d<MesonOperatorInsertion_t> inser;  // outside array is over gammas
  };


  //! Serialize generalized operator object
  multi1d<Complex> MesonOperator_t::serialize()
  {
    int inser_size = inser.size();
    int op_size2   = inser[0].op.size2();
    int op_size1   = inser[0].op.size1();
    int elem_size2 = inser[0].op(0,0).elem.size2();
    int elem_size1 = inser[0].op(0,0).elem.size1();

    // dreadful hack - use a complex to hold an int
    Complex inser_sizes, op_sizes, elem_sizes;
    inser_sizes = cmplx(Real(inser.size()), Real(zero));
    op_sizes    = cmplx(Real(op_size2), Real(op_size1));
    elem_sizes  = cmplx(Real(elem_size2), Real(elem_size1));

    multi1d<Complex> mesprop_1d(3 + inser_size*op_size2*op_size1*elem_size2*elem_size1);

    int cnt = 0;

    mesprop_1d[cnt++] = inser_sizes;
    mesprop_1d[cnt++] = op_sizes;
    mesprop_1d[cnt++] = elem_sizes;

    for(int g=0; g < inser.size(); ++g)             // inser
      for(int i=0; i < inser[g].op.size2(); ++i)             // op_l
	for(int j=0; j < inser[g].op.size1(); ++j)           // op_r
	  for(int k=0; k < inser[g].op(i,j).elem.size2(); ++k)    // elem_l
	    for(int l=0; l < inser[g].op(i,j).elem.size1(); ++l)  // elem_r
	      mesprop_1d[cnt++] = inser[g].op(i,j).elem(k,l);

    return mesprop_1d;
  }


  //! MesonOperator header writer
  void write(XMLWriter& xml, const string& path, const MesonOperator_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "mom2_max", param.mom2_max);
    write(xml, "j_decay", param.j_decay);
    write(xml, "seed_l", param.seed_l);
    write(xml, "seed_r", param.seed_r);
    write(xml, "smearing_l", param.smearing_l);
    write(xml, "smearing_r", param.smearing_r);

    if( path != "." )
      pop(xml);
  }



  //--------------------------------------------------------------
  // Function call
  void 
  InlineStochMeson::operator()(unsigned long update_no,
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
  InlineStochMeson::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

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
      QDPIO::cerr << InlineStochMesonEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineStochMesonEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "stoch_meson");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineStochMesonEnv::name << ": Stochastic Meson Operator" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

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
    snoop.start();
    multi1d<QuarkSourceSolutions_t>  quarks(params.named_obj.prop.op.size());
    QDPIO::cout << "num_quarks= " << params.named_obj.prop.op.size() << endl;

    try
    {
      QDPIO::cout << "quarks.size= " << quarks.size() << endl;
      for(int n=0; n < quarks.size(); ++n)
      {
	QDPIO::cout << "Attempt to read solutions for source number=" << n << endl;
	quarks[n].dilutions.resize(params.named_obj.prop.op[n].soln_files.size());

	QDPIO::cout << "dilutions.size= " << quarks[n].dilutions.size() << endl;
	for(int i=0; i < quarks[n].dilutions.size(); ++i)
	{
	  XMLReader file_xml, record_xml;

	  QDPIO::cout << "reading file= " << params.named_obj.prop.op[n].soln_files[i] << endl;
	  QDPFileReader from(file_xml, params.named_obj.prop.op[n].soln_files[i], QDPIO_SERIAL);
	  read(from, record_xml, quarks[n].dilutions[i].soln);
	  close(from);
	
	  read(record_xml, "/Propagator/PropSource", quarks[n].dilutions[i].source_header);
	  read(record_xml, "/Propagator/ForwardProp", quarks[n].dilutions[i].prop_header);
	}
      }
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting headers: " << e << endl;
      QDP_abort(1);
    }
    snoop.stop();

    QDPIO::cout << "Sources and solutions successfully read: time= "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;



    //
    // Check for each quark source that the solutions have their diluted
    // on every site only once
    //
    snoop.start();

    try
    {
      push(xml_out, "Norms");
      for(int n=0; n < quarks.size(); ++n)
      {
	bool first = true;
	int  N;
	LatticeFermion quark_noise;      // noisy source on entire lattice

	for(int i=0; i < quarks[n].dilutions.size(); ++i)
	{
	  std::istringstream  xml_s(quarks[n].dilutions[i].source_header.source.xml);
	  XMLReader  sourcetop(xml_s);
//	QDPIO::cout << "Source = " << quarks[n].dilutions[i].source_header.source.id << endl;

	  if (quarks[n].dilutions[i].source_header.source.id != DiluteZNQuarkSourceConstEnv::getName())
	  {
	    QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::getName() << endl;
	    QDP_abort(1);
	  }

	  QDPIO::cout << "Quark num= " << n << "  dilution num= " << i << endl;

	  // Manually create the params so I can peek into them and use the source constructor
	  DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
							 quarks[n].dilutions[i].source_header.source.path);
	  DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
	  if (first) 
	  {
	    first = false;

	    quarks[0].j_decay = srcParams.j_decay;

	    // Grab N
	    N = srcParams.N;

	    // Set the seed to desired value
	    quarks[n].seed = srcParams.ran_seed;
	    QDP::RNG::setrn(quarks[n].seed);

	    // Create the noisy quark source on the entire lattice
	    zN_src(quark_noise, N);
	  }

	  // The seeds must always agree - here the seed is the unique id of the source
	  if ( toBool(srcParams.ran_seed != quarks[n].seed) )
	  {
	    QDPIO::cerr << "quark source=" << n << "  dilution=" << i << " seed does not match" << endl;
	    QDP_abort(1);
	  }

	  // The N's must always agree
	  if ( toBool(srcParams.N != N) )
	  {
	    QDPIO::cerr << "quark source=" << n << "  dilution=" << i << " N does not match" << endl;
	    QDP_abort(1);
	  }

#if 0
	  // Use a trick here, create the source and subtract it from the global noisy
	  // Check at the end that the global noisy is zero everywhere.
	  // NOTE: the seed will be set every call
	  quarks[n].dilutions[i].source = srcConst(u);
	  quark_noise -= quarks[n].dilutions[i].source;

	  // Diagnostic
	  {
	    // Keep a copy of the phases with NO momenta

	    multi1d<Double> source_corr = sumMulti(localNorm2(quarks[n].dilutions[i].source), 
						   phases_nomom.getSet());

	    multi1d<Double> soln_corr = sumMulti(localNorm2(quarks[n].dilutions[i].soln), 
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

#if 0
	Double dcnt = norm2(quark_noise, phases );
	if (toDouble(dcnt) != 0.0)  // problematic - seems to work with unnormalized sources 
	{
	  QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	  QDP_abort(1);
	}

#endif
      } // end for n

      pop(xml_out);
    } // end try
    catch(const std::string& e) 
    {
      QDPIO::cerr << ": Caught Exception creating source: " << e << endl;
      QDP_abort(1);
    }

    snoop.stop();

    QDPIO::cout << "Sources saturated: time= "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;


    //
    // Meson operators
    //
    int j_decay = quarks[0].j_decay;

    // Initialize the slow Fourier transform phases
    SftMom phases(params.param.mom2_max, false, j_decay);
    
    // Length of lattice in decay direction
    int length = phases.numSubsets();


    if (quarks.size() != 2)
    {
      QDPIO::cerr << "expecting 2 quarks but have num quarks= " << quarks.size() << endl;
      QDP_abort(1);
    }


    // Operator A
    snoop.start();
    MesonOperator_t  meson_opA;
    meson_opA.mom2_max    = params.param.mom2_max;
    meson_opA.j_decay     = j_decay;
    meson_opA.seed_l      = quarks[1].seed;
    meson_opA.seed_r      = quarks[0].seed;
    meson_opA.smearing_l  = params.source_smearing.source.xml;
    meson_opA.smearing_r  = params.source_smearing.source.xml;
    meson_opA.inser.resize(Ns*Ns);

    // Sanity check
    if ( toBool(meson_opA.seed_l == meson_opA.seed_r) )
    {
      QDPIO::cerr << "meson op seeds are the same" << endl;
      QDP_abort(1);
    }

    // Construct operator A
    try
    {
      int G5 = Ns*Ns-1;

      std::istringstream  xml_s(params.source_smearing.source.xml);
      XMLReader  sourcetop(xml_s);
      QDPIO::cout << "Source = " << params.source_smearing.source.id << endl;

      Handle< QuarkSourceSink<LatticeFermion> >
	sourceSmearing(TheFermSourceSmearingFactory::Instance().createObject(
			 params.source_smearing.source.id,
			 sourcetop,
			 params.source_smearing.source.path,
			 u));

      push(xml_out, "OperatorA");

      // Source smear all the sources up front
      multi1d<LatticeFermion> smeared_sources(quarks[1].dilutions.size());
      for(int i=0; i < smeared_sources.size(); ++i)
      {
	smeared_sources[i] = quarks[1].dilutions[i].source;
	(*sourceSmearing)(smeared_sources[i]);
      }

      // Source smear all the solutions up front
      multi1d<LatticeFermion> smeared_solns(quarks[0].dilutions.size());
      for(int j=0; j < smeared_solns.size(); ++j)
      {
	smeared_solns[j] = quarks[0].dilutions[j].soln;
	(*sourceSmearing)(smeared_solns[j]);
      }

      QDPIO::cout << "source smearings done" << endl;

      for(int gamma_value=0; gamma_value < Ns*Ns; ++gamma_value)
      {
	meson_opA.inser[gamma_value].op.resize(smeared_sources.size(), smeared_solns.size());

	for(int i=0; i < smeared_sources.size(); ++i)
	{
	  for(int j=0; j < smeared_solns.size(); ++j)
	  {
	    // Optimize by restricting operations to source time slice
	    LatticeComplex corr_fn = zero;
	    corr_fn[phases.getSet()[quarks[1].dilutions[i].source_header.t_source]] = 
	      localInnerProduct(smeared_sources[i], Gamma(gamma_value) * smeared_solns[j]);
	    meson_opA.inser[gamma_value].op(i,j).elem = phases.sft(corr_fn);
	  } // end for j
	} // end for i
      } // end for g

      pop(xml_out); // OperatorA
    } // opA
    catch(const std::string& e) 
    {
      QDPIO::cerr << ": Caught Exception creating source smearing: " << e << endl;
      QDP_abort(1);
    }
    catch(...)
    {
      QDPIO::cerr << ": Caught generic exception creating source smearing" << endl;
      QDP_abort(1);
    }

    snoop.stop();

    QDPIO::cout << "Operator A computed: time= "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;


    // Operator B
    snoop.start();
    MesonOperator_t  meson_opB;
    meson_opB.mom2_max    = params.param.mom2_max;
    meson_opB.j_decay     = j_decay;
    meson_opB.seed_l      = quarks[0].seed;
    meson_opB.seed_r      = quarks[1].seed;
    meson_opB.smearing_l  = params.sink_smearing.sink.xml;
    meson_opB.smearing_r  = params.sink_smearing.sink.xml;
    meson_opB.inser.resize(Ns*Ns);

    // Sanity check
    if ( toBool(meson_opB.seed_l == meson_opB.seed_r) )
    {
      QDPIO::cerr << "meson op seeds are the same" << endl;
      QDP_abort(1);
    }

    // Construct operator B
    {
      int G5 = Ns*Ns-1;

      std::istringstream  xml_s(params.sink_smearing.sink.xml);
      XMLReader  sinktop(xml_s);
      QDPIO::cout << "Sink = " << params.sink_smearing.sink.id << endl;

      Handle< QuarkSourceSink<LatticeFermion> >
	sinkSmearing(TheFermSinkSmearingFactory::Instance().createObject(
		       params.sink_smearing.sink.id,
		       sinktop,
		       params.sink_smearing.sink.path,
		       u));

      push(xml_out, "OperatorB");

      // Sink smear all the sources up front
      multi1d<LatticeFermion> smeared_sources(quarks[0].dilutions.size());
      for(int j=0; j < smeared_sources.size(); ++j)
      {
	smeared_sources[j] = quarks[0].dilutions[j].source;
	(*sinkSmearing)(smeared_sources[j]);
      }

      // Sink smear all the solutions up front
      multi1d<LatticeFermion> smeared_solns(quarks[1].dilutions.size());
      for(int i=0; i < smeared_solns.size(); ++i)
      {
	smeared_solns[i] = quarks[1].dilutions[i].soln;
	(*sinkSmearing)(smeared_solns[i]);
      }

      QDPIO::cout << "sink smearings done" << endl;

      for(int gamma_value=0; gamma_value < Ns*Ns; ++gamma_value)
      {
	meson_opB.inser[gamma_value].op.resize(smeared_sources.size(), smeared_solns.size());

	for(int j=0; j < smeared_sources.size(); ++j)
	{
	  for(int i=0; i < smeared_solns.size(); ++i)
	  {
	    // Optimize by restricting operations to source time slice
	    LatticeComplex corr_fn = zero;
	    corr_fn[phases.getSet()[quarks[0].dilutions[j].source_header.t_source]] = 
	      localInnerProduct(smeared_sources[j], Gamma(gamma_value) * smeared_solns[i]);
	    meson_opB.inser[gamma_value].op(j,i).elem = phases.sft(corr_fn);
	  } // end for i
	} // end for j
      } // end for g

      pop(xml_out); // OperatorB
    } // opB

    snoop.stop();

    QDPIO::cout << "Operator B computed: time= "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;


    // Save the operators
    // ONLY SciDAC output format is supported!
    snoop.start();
    {
      XMLBufferWriter file_xml;
      push(file_xml, "meson_operator");
      write(file_xml, "Config_info", gauge_xml);
      pop(file_xml);

      QDPFileWriter to(file_xml, params.named_obj.prop.op_file,     // are there one or two files???
		       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

      // Write the scalar data
      {
	XMLBufferWriter record_xml;
	write(record_xml, "SourceMesonOperator", meson_opA);
	write(to, record_xml, meson_opA.serialize());
      }

      // Write the scalar data
      {
	XMLBufferWriter record_xml;
	write(record_xml, "SinkMesonOperator", meson_opB);
	write(to, record_xml, meson_opB.serialize());
      }

      close(to);
    }

    snoop.stop();

    QDPIO::cout << "Operators written: time= "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    // Close the namelist output file XMLDAT
    pop(xml_out);     // StochMeson

    snoop.stop();
    QDPIO::cout << InlineStochMesonEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineStochMesonEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

}
