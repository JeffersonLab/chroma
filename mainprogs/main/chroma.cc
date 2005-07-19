// $Id: chroma.cc,v 1.5 2005-07-19 22:27:09 edwards Exp $
/*! \file
 *  \brief Main program to run all measurement codes.
 */

#include "chroma.h"

using namespace Chroma;
extern "C" { 
 void _mcleanup();
};

/*
 * Input 
 */
struct Params_t
{
  multi1d<int>    nrow;
  std::string     inline_measurement_xml;
};

struct Inline_input_t
{
  Params_t        param;
  Cfg_t           cfg;
  QDP::Seed       rng_seed;
};


void read(XMLReader& xml, const std::string& path, Params_t& p) 
{
  XMLReader paramtop(xml, path);
  read(paramtop, "nrow", p.nrow);

  XMLReader measurements_xml(paramtop, "InlineMeasurements");
  std::ostringstream inline_os;
  measurements_xml.print(inline_os);
  p.inline_measurement_xml = inline_os.str();
  QDPIO::cout << "InlineMeasurements are: " << endl;
  QDPIO::cout << p.inline_measurement_xml << endl;
}


void read(XMLReader& xml, const std::string& path, Inline_input_t& p) 
{
  try {
    XMLReader paramtop(xml, path);
      
    read(paramtop, "Param", p.param);
    read(paramtop, "Cfg", p.cfg);

    if (paramtop.count("RNG") > 0)
      read(paramtop, "RNG", p.rng_seed);
    else
      p.rng_seed = 11;     // default seed
  }
  catch( const std::string& e ) 
  {
    QDPIO::cerr << "Error reading XML : " << e << endl;
    QDP_abort(1);
  }
}

 

bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= InlineAggregateEnv::registered;

  return foo;
}

//! Main program to run all measurement codes
/*! \defgroup chromamain Main program to run all measurement codes.
 *  \ingroup main
 *
 * Main program to run all measurement codes.
 */

int main(int argc, char *argv[]) 
{
  // Chroma Init stuff
  Chroma::initialize(&argc, &argv);
  
  linkageHack();

  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Input parameter structure
  Inline_input_t  input;
  read(xml_in, "/chroma", input);

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "chroma");

  // Write out the input
  write(xml_out, "Input", xml_in);

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  proginfo(xml_out);    // Print out basic program info

  // Initialise the RNG
  QDP::RNG::setrn(input.rng_seed);
  write(xml_out,"RNG", input.rng_seed);

  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Start up the gauge field
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  XMLBufferWriter config_xml;
  config_xml << gauge_xml;

  // Write out the config header
  push(xml_out, "Config_info");
  write(xml_out, "file_xml", gauge_file_xml);
  write(xml_out, "gauge_xml", gauge_xml);
  pop(xml_out);
  
  // Get the measurements
  std::istringstream Measurements_is(input.param.inline_measurement_xml);
  XMLReader MeasXML(Measurements_is);
  multi1d < Handle< AbsInlineMeasurement > > the_measurements;

  try {
    read(MeasXML, "/InlineMeasurements", the_measurements);
  }
  catch(const std::string& e) {
    QDPIO::cerr << "Caught Exception: " << e << endl;
    QDP_abort(1);
  }

  QDPIO::cout << "There are " << the_measurements.size() << " measurements " << endl;

  // Measure inline observables 
  push(xml_out, "InlineObservables");

  QDPIO::cout << "Doing " << the_measurements.size() 
	      <<" measurements" << endl;
  unsigned long cur_update = 0;
  for(int m=0; m < the_measurements.size(); m++) 
  {
    AbsInlineMeasurement& the_meas = *(the_measurements[m]);
    if( cur_update % the_meas.getFrequency() == 0 ) 
    {
      // Caller writes elem rule
      push(xml_out, "elem");
      the_meas(u, config_xml, cur_update, xml_out);
      pop(xml_out); 
    }
  }
  pop(xml_out); // pop("InlineObservables");
  pop(xml_out);

  Chroma::finalize();

  
  exit(0);
}

