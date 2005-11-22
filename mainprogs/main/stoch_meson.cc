// $Id: stoch_meson.cc,v 1.5 2005-11-22 22:02:53 edwards Exp $
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
    std::string          op_file;
  };

  multi1d<Operator_t>  op;
};


//! Mega-structure of all input
struct StochMeson_input_t
{
  Param_t       param;
  PropSource_t  source_smearing;
  PropSink_t    sink_smearing;
  Prop_t        prop;
  Cfg_t         cfg;

  multi1d<int> nrow;
};


//! Operator parameters
void read(XMLReader& xml, const string& path, Prop_t::Operator_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "soln_files", input.soln_files);
  read(inputtop, "operator_file", input.op_file);
}


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

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
    LatticeFermion  source;
    LatticeFermion  soln;
    PropSource_t    source_header;
    ChromaProp_t    prop_header;
  };

  int   j_decay;
  Seed  seed;
  multi1d<QuarkSolution_t>  dilutions;
};


//! Meson operator
struct MesonOperator_t
{
  //! Meson operator
  struct MesonOperatorElement_t
  {
    multi2d<DComplex> elem;
  };

  int   j_decay;
  Seed  seed_l;
  Seed  seed_r;
  multi2d<MesonOperatorElement_t> op;
};


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

  QDPIO::cout << " STOCHMESON: stochastically estimate a meson operator" << endl;

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
  QDPIO::cout << "Sources and solutions successfully read\n" << endl;


  //
  // Check for each quark source that the solutions have their diluted
  // on every site only once
  //
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
  MesonOperator_t  meson_opA;
  meson_opA.j_decay = j_decay;
  meson_opA.seed_l  = quarks[1].seed;
  meson_opA.seed_r  = quarks[0].seed;
  meson_opA.op.resize(quarks[1].dilutions.size(), quarks[0].dilutions.size());

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
    int gamma = 0;   // need to understand this convention - I thought I should use G5 for a pion

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

    // Source smear all the solutions up front
    multi1d<LatticeFermion> smeared_solns(meson_opA.op.size1());
    for(int j=0; j < meson_opA.op.size1(); ++j)
    {
      smeared_solns[j] = quarks[0].dilutions[j].soln;
      (*sourceSmearing)(smeared_solns[j]);
    }

    // Could have some optimizations for time slice dilutions
    for(int i=0; i < meson_opA.op.size2(); ++i)
    {
      LatticeFermion source_tmp = quarks[1].dilutions[i].source;
      (*sourceSmearing)(source_tmp);

      for(int j=0; j < meson_opA.op.size1(); ++j)
      {
	LatticeFermion z = Gamma(G5) * (Gamma(gamma) * smeared_solns[j]); // do lots of stuff here
	LatticeComplex corr_fn = localInnerProduct(source_tmp, z);
	meson_opA.op(i,j).elem = phases.sft(corr_fn);

#if 0
	// Diagnostic
	{
	  push(xml_out, "elem");
	  write(xml_out, "i", i);
	  write(xml_out, "j", j);
	  write(xml_out, "opA", meson_opA.op(i,j).elem[0]);
	  pop(xml_out);
	}
#endif
      } // end for j
    } // end for i

    pop(xml_out); // OperatorA
  } // opA
  catch(const std::string& e) 
  {
    QDPIO::cerr << ": Caught Exception creating source smearing: " << e << endl;
    QDP_abort(1);
  }


  // Operator B
  MesonOperator_t  meson_opB;
  meson_opB.j_decay = j_decay;
  meson_opB.seed_l  = quarks[0].seed;
  meson_opB.seed_r  = quarks[1].seed;
  meson_opB.op.resize(quarks[0].dilutions.size(), quarks[1].dilutions.size());

  // Sanity check
  if ( toBool(meson_opB.seed_l == meson_opB.seed_r) )
  {
    QDPIO::cerr << "meson op seeds are the same" << endl;
    QDP_abort(1);
  }

  // Construct operator B
  {
    int G5 = Ns*Ns-1;
    int gamma = 0;   // need to understand this convention - I thought I should use G5 for a pion

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

    // Sink smear all the solutions up front
    multi1d<LatticeFermion> smeared_solns(meson_opB.op.size1());
    for(int i=0; i < meson_opB.op.size1(); ++i)
    {
      smeared_solns[i] = quarks[1].dilutions[i].soln;
      (*sinkSmearing)(smeared_solns[i]);
    }

    // Could have some optimizations for time slice dilutions
    for(int j=0; j < meson_opB.op.size2(); ++j)
    {
      LatticeFermion source_tmp = quarks[0].dilutions[j].source;
      (*sinkSmearing)(source_tmp);

      for(int i=0; i < meson_opB.op.size1(); ++i)
      {
	LatticeFermion z = Gamma(G5) * (Gamma(gamma) * smeared_solns[i]); // do lots of stuff here
	LatticeComplex corr_fn = localInnerProduct(source_tmp, z);
	meson_opB.op(j,i).elem = phases.sft(corr_fn);

#if 0
	// Diagnostic
	{
	  push(xml_out, "elem");
	  write(xml_out, "j", j);
	  write(xml_out, "i", i);
	  write(xml_out, "opB", meson_opB.op(j,i).elem[0]);
	  pop(xml_out);
	}
#endif
      } // end for i
    } // end for j

    pop(xml_out); // OperatorB
  } // opB



  //
  // For now, compute a pion correlator. This part should be moved to ADAT
  //
  {
    push(xml_out, "Mesons");

    // Loop over sink momenta
    XMLArrayWriter xml_sink_mom(xml_out,phases.numMom());
    push(xml_sink_mom, "momenta");

    // The momenta loop is not handled correctly yet
    for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
    {
      push(xml_sink_mom);
      write(xml_sink_mom, "sink_mom_num", sink_mom_num);
      write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num));

      multi1d<Double> mesprop(length);
      multi1d<DComplex> cmesprop(length);
      mesprop = zero;
      cmesprop = zero;

      for(int dt=0; dt < length; ++dt) 
	for(int t=0; t < length; ++t) 
	  for(int i=0; i < meson_opA.op.size2(); ++i)
	    for(int j=0; j < meson_opA.op.size1(); ++j)
	    {
	      cmesprop[dt] += (meson_opA.op(i,j).elem[sink_mom_num][t] *
			       meson_opB.op(j,i).elem[sink_mom_num][(t+dt) % length]);

	      mesprop[dt] += real(meson_opA.op(i,j).elem[sink_mom_num][t] *
				  meson_opB.op(j,i).elem[sink_mom_num][(t+dt) % length]);
	    }

      write(xml_sink_mom, "cmesprop", cmesprop);
      write(xml_sink_mom, "mesprop", mesprop);
      pop(xml_sink_mom);
    } // end for(sink_mom_num)
      
    pop(xml_sink_mom);
    pop(xml_out); // Mesons
  }


  // Close the namelist output file XMLDAT
  pop(xml_out);     // StochMeson

  xml_in.close();
  xml_out.close();

  END_CODE();

// Time to bolt
  QDP_finalize();

  exit(0);
}
