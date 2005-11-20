// $Id: stoch_meson.cc,v 1.2 2005-11-20 18:28:07 edwards Exp $
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
  multi1d<int> nrow;
};

//! Propagators
struct Prop_t
{
  multi1d< multi1d<std::string> >  soln_files;
};


//! Mega-structure of all input
struct StochMeson_input_t
{
  Param_t     param;
  Prop_t      prop;
  Cfg_t       cfg;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "soln_files", input.soln_files);
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

  read(paramtop, "nrow", param.nrow);
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

    // Read in the propagator(s) info
    read(inputtop, "Prop", input.prop);

    // Read the gauge field
    read(inputtop, "Cfg", input.cfg);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading stoch_meson data: " << e << endl;
    QDP_abort(1);
  }
}

//--------------------------------------------------------------


//! Structure holding solutions
struct QuarkSolution_t
{
  LatticeFermion  quark;
  PropSource_t    source_header;
  ChromaProp_t    prop_header;
};

//! Structure holding a source and its solutions
struct QuarkSourceSolutions_t
{
  multi1d<QuarkSolution_t>  dilutions;
};



//! Stochastically estimate meson operators
/*! \defgroup stochastic Stochastically estimate meson operators
 *  \ingroup main
 *
 * Main program to stochastically estimate meson operators
 */

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

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
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " STOCHMESON: stochastically estimate a meson operator" << endl;

  // Instantiate XML writer for XMLDAT
  // XMLFileWriter  xml_out(Chroma::getXMLOutputFileName());
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "stoch_meson");

  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Start up the gauge field
  QDPIO::cout << "Attempt to read gauge field" << endl;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);
  QDPIO::cout << "Gauge field successfully read" << endl;

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();

  // Save current seed
  Seed ran_seed;
  QDP::RNG::savern(ran_seed);

  //
  // Read the source and solutions
  //
  multi1d<QuarkSourceSolutions_t>  quarks(input.prop.soln_files.size());
  QDPIO::cout << "num_quarks= " << input.prop.soln_files.size() << endl;

  try
  {
    QDPIO::cout << "quarks.size= " << quarks.size() << endl;
    for(int n=0; n < quarks.size(); ++n)
    {
      QDPIO::cout << "Attempt to read solutions for source number=" << n << endl;
      quarks[n].dilutions.resize(input.prop.soln_files[n].size());

      QDPIO::cout << "dilutions.size= " << quarks[n].dilutions.size() << endl;
      for(int i=0; i < quarks[n].dilutions.size(); ++i)
      {
	XMLReader file_xml, record_xml;

	QDPIO::cout << "reading file= " << input.prop.soln_files[n][i] << endl;
	QDPFileReader from(file_xml, input.prop.soln_files[n][i], QDPIO_SERIAL);
	read(from, record_xml, quarks[n].dilutions[i].quark);
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


#if 0
  // Sanity check - write out the norm2 of the forward prop in the Nd-1 direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);

    multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator1), 
						 phases.getSet());

    push(xml_out, "Forward_prop_correlator1");
    write(xml_out, "forward_prop_corr", forward_prop_corr);
    pop(xml_out);
  }

  // Save prop input
  push(xml_out, "Propagators");
  write(xml_out, "ForwardProp", soln_prop_headers);
  write(xml_out, "PropSource", soln_source_headers);
  pop(xml_out);
#endif


  //
  // Check for each quark source that the solutions have their diluted
  // on every site only once
  //
  try
  {
    for(int n=0; n < quarks.size(); ++n)
    {
      bool first = true;
      int  N;
      Seed source_seed;
      LatticeFermion quark_noise;      // noisy source on entire lattice

      for(int i=0; i < quarks[n].dilutions.size(); ++i)
      {
	std::istringstream  xml_s(quarks[n].dilutions[i].source_header.source);
	XMLReader  sourcetop(xml_s);
	const string source_path = "/Source";
	QDPIO::cout << "Source = " << quarks[n].dilutions[i].source_header.source_type << endl;

	if (quarks[n].dilutions[i].source_header.source_type != DiluteZNQuarkSourceConstEnv::name)
	{
	  QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
	  QDP_abort(1);
	}

	// Manually create the params so I can peek into them and use the source constructor
	DiluteZNQuarkSourceConstEnv::Params       srcParams(sourcetop, source_path);
	DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
	if (first) 
	{
	  first = false;

	  // Grab N
	  N = srcParams.N;

	  // Set the seed to desired value
	  source_seed = srcParams.ran_seed;
	  QDP::RNG::setrn(source_seed);

	  // Create the noisy quark source on the entire lattice
	  zN_src(quark_noise, N);
	}

	// The seeds must always agree - here the seed is the unique id of the source
	if ( toBool(srcParams.ran_seed != source_seed) )
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
	quark_noise -= srcConst(u);

      } // end for i

      Double dcnt = norm2(quark_noise);
      if (toDouble(dcnt) != 0.0)  // problematic - seems to work with unnormalized sources 
      {
	QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	QDP_abort(1);
      }

    } // end for n
  } // end try
  catch(const std::string& e) 
  {
    QDPIO::cerr << ": Caught Exception creating source: " << e << endl;
    QDP_abort(1);
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
