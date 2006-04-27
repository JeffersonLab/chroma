// $Id: stoch_meson.cc,v 3.1 2006-04-27 02:36:00 edwards Exp $
/*! \file
 * \brief Stochastically estimate a meson operator
 *
 * This is the first crude attempt at constructing a meson operator
 * from stochastic sources. In this case, the sources and solutions
 * must come from a "diluted" source.
 */

#include "chroma.h"

using namespace Chroma;


/*
 * Input 
 */
//! Parameters for running program
struct Param_t
{
  int mom2_max;            // (mom)^2 <= mom2_max
};


//! Propagators
struct Prop_t
{
  //! Operators
  struct Operator_t
  {
    multi1d<std::string> soln_files;
  };

  std::string          op_file;
  multi1d<Operator_t>  op;
};


//! Mega-structure of all input
struct StochMeson_input_t
{
  Param_t            param;
  PropSourceSmear_t  source_smearing;
  PropSinkSmear_t    sink_smearing;
  Prop_t             prop;
  Cfg_t              cfg;

  multi1d<int> nrow;
};


//! Operator parameters
void read(XMLReader& xml, const string& path, Prop_t::Operator_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "soln_files", input.soln_files);
}


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "operator_file", input.op_file);
  read(inputtop, "operator", input.op);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Param_t& param)
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
void read(XMLReader& xml, const string& path, StochMeson_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Source smearing
    read(inputtop, "SourceSmearing", input.source_smearing);

    // Sink smearing
    read(inputtop, "SinkSmearing", input.sink_smearing);

    // Read in the propagator(s) info
    read(inputtop, "Prop", input.prop);

    // Read the gauge field
    read(inputtop, "Cfg", input.cfg);

    // Lattice size
    read(inputtop, "nrow", input.nrow);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading stoch_meson data: " << e << endl;
    QDP_abort(1);
  }
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




bool linkageHack(void)
{
  bool foo = true;

  // Source and sink smearing
  foo &= QuarkSourceSmearingEnv::registered;
  foo &= QuarkSinkSmearingEnv::registered;

  return foo;
}

//! Stochastically estimate meson operators
/*! \defgroup stochastic Stochastically estimate meson operators
 *  \ingroup main
 *
 * Main program to stochastically estimate meson operators
 */

int main(int argc, char *argv[])
{
  START_CODE();

  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  linkageHack();

  XMLReader xml_in;

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  // Input parameter structure
  StochMeson_input_t  input;

  try
  {
    // Read input file
    xml_in.open(Chroma::getXMLInputFileName());
    read(xml_in, "/stoch_meson", input);
  }
  catch(const std::string& e) 
  {
    QDPIO::cerr << "STOCH_MESON: Caught Exception reading XML: " << e << endl;
    QDP_abort(1);
  }
  catch(...)
  {
    QDPIO::cerr << "STOCH_MESON: caught generic exception reading XML" << endl;
    QDP_abort(1);
  }

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.nrow);
  Layout::create();

  QDPIO::cout << " STOCH_MESON: stochastically estimate a meson operator" << endl;

  // Instantiate XML writer for XMLDAT
  // XMLFileWriter  xml_out(Chroma::getXMLOutputFileName());
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "stoch_meson");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();

  // Start up the config
  StopWatch swatch;
  swatch.reset();
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Start up the gauge field
  QDPIO::cout << "Attempt to read gauge field" << endl;
  swatch.start();
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);
  swatch.stop();

  QDPIO::cout << "Gauge field successfully read: time= " 
	      << swatch.getTimeInSeconds() 
	      << " secs" << endl;

  // Write out the config header
  push(xml_out, "Config_info");
  write(xml_out, "file_xml", gauge_file_xml);
  write(xml_out, "gauge_xml", gauge_xml);
  pop(xml_out);
  
  // Calculate some gauge invariant observables
  MesPlq(xml_out, "Observables", u);

  xml_out.flush();

  // Save current seed
  Seed ran_seed;
  QDP::RNG::savern(ran_seed);

  //
  // Read the source and solutions
  //
  swatch.start();
  multi1d<QuarkSourceSolutions_t>  quarks(input.prop.op.size());
  QDPIO::cout << "num_quarks= " << input.prop.op.size() << endl;

  try
  {
    QDPIO::cout << "quarks.size= " << quarks.size() << endl;
    for(int n=0; n < quarks.size(); ++n)
    {
      QDPIO::cout << "Attempt to read solutions for source number=" << n << endl;
      quarks[n].dilutions.resize(input.prop.op[n].soln_files.size());

      QDPIO::cout << "dilutions.size= " << quarks[n].dilutions.size() << endl;
      for(int i=0; i < quarks[n].dilutions.size(); ++i)
      {
	XMLReader file_xml, record_xml;

	QDPIO::cout << "reading file= " << input.prop.op[n].soln_files[i] << endl;
	QDPFileReader from(file_xml, input.prop.op[n].soln_files[i], QDPIO_SERIAL);
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
  swatch.stop();

  QDPIO::cout << "Sources and solutions successfully read: time= "
	      << swatch.getTimeInSeconds() 
	      << " secs" << endl;



  //
  // Check for each quark source that the solutions have their diluted
  // on every site only once
  //
  swatch.start();

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
	std::istringstream  xml_s(quarks[n].dilutions[i].source_header.source);
	XMLReader  sourcetop(xml_s);
	const string source_path = "/Source";
//	QDPIO::cout << "Source = " << quarks[n].dilutions[i].source_header.source_type << endl;

	if (quarks[n].dilutions[i].source_header.source_type != DiluteZNQuarkSourceConstEnv::name)
	{
	  QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
	  QDP_abort(1);
	}

	QDPIO::cout << "Quark num= " << n << "  dilution num= " << i << endl;

	// Manually create the params so I can peek into them and use the source constructor
	DiluteZNQuarkSourceConstEnv::Params       srcParams(sourcetop, source_path);
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

	// Use a trick here, create the source and subtract it from the global noisy
	// Check at the end that the global noisy is zero everywhere.
	// NOTE: the seed will be set every call
	quarks[n].dilutions[i].source = srcConst(u);
	quark_noise -= quarks[n].dilutions[i].source;

#if 0
	// Diagnostic
	{
	  // Keep a copy of the phases with NO momenta
	  SftMom phases_nomom(0, true, quarks[n].dilutions[i].source_header.j_decay);

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

      Double dcnt = norm2(quark_noise);
      if (toDouble(dcnt) != 0.0)  // problematic - seems to work with unnormalized sources 
      {
	QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	QDP_abort(1);
      }

    } // end for n

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
  int j_decay = quarks[0].j_decay;

  // Initialize the slow Fourier transform phases
  SftMom phases(input.param.mom2_max, false, j_decay);
    
  // Length of lattice in decay direction
  int length = phases.numSubsets();


  if (quarks.size() != 2)
  {
    QDPIO::cerr << "expecting 2 quarks but have num quarks= " << quarks.size() << endl;
    QDP_abort(1);
  }


  // Operator A
  swatch.start();
  MesonOperator_t  meson_opA;
  meson_opA.mom2_max    = input.param.mom2_max;
  meson_opA.j_decay     = j_decay;
  meson_opA.seed_l      = quarks[1].seed;
  meson_opA.seed_r      = quarks[0].seed;
  meson_opA.smearing_l  = input.source_smearing.source;
  meson_opA.smearing_r  = input.source_smearing.source;
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

    std::istringstream  xml_s(input.source_smearing.source);
    XMLReader  sourcetop(xml_s);
    const string source_path = "/Source";
    QDPIO::cout << "Source = " << input.source_smearing.source_type << endl;

    Handle< QuarkSourceSink<LatticeFermion> >
      sourceSmearing(TheFermSourceSmearingFactory::Instance().createObject(
		       input.source_smearing.source_type,
		       sourcetop,
		       source_path,
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

  swatch.stop();

  QDPIO::cout << "Operator A computed: time= "
	      << swatch.getTimeInSeconds() 
	      << " secs" << endl;


  // Operator B
  swatch.start();
  MesonOperator_t  meson_opB;
  meson_opB.mom2_max    = input.param.mom2_max;
  meson_opB.j_decay     = j_decay;
  meson_opB.seed_l      = quarks[0].seed;
  meson_opB.seed_r      = quarks[1].seed;
  meson_opB.smearing_l  = input.sink_smearing.sink;
  meson_opB.smearing_r  = input.sink_smearing.sink;
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

    std::istringstream  xml_s(input.sink_smearing.sink);
    XMLReader  sinktop(xml_s);
    const string sink_path = "/Sink";
    QDPIO::cout << "Sink = " << input.sink_smearing.sink_type << endl;

    Handle< QuarkSourceSink<LatticeFermion> >
      sinkSmearing(TheFermSinkSmearingFactory::Instance().createObject(
		     input.sink_smearing.sink_type,
		     sinktop,
		     sink_path,
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

  swatch.stop();

  QDPIO::cout << "Operator B computed: time= "
	      << swatch.getTimeInSeconds() 
	      << " secs" << endl;


  // Save the operators
  // ONLY SciDAC output format is supported!
  swatch.start();
  {
    XMLBufferWriter file_xml;
    push(file_xml, "meson_operator");
    write(file_xml, "Config_info", gauge_xml);
    pop(file_xml);

    QDPFileWriter to(file_xml, input.prop.op_file,     // are there one or two files???
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

  swatch.stop();

  QDPIO::cout << "Operators written: time= "
	      << swatch.getTimeInSeconds() 
	      << " secs" << endl;

  // Close the namelist output file XMLDAT
  pop(xml_out);     // StochMeson

  xml_in.close();
  xml_out.close();

  snoop.stop();
  QDPIO::cout << "stoch_meson: total time= "
	      << snoop.getTimeInSeconds() 
	      << " secs" << endl;

  END_CODE();

// Time to bolt
  QDP_finalize();

  exit(0);
}
