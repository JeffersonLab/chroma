// $Id: t_inline.cc,v 1.1 2005-02-08 03:17:22 edwards Exp $
// Test driver for inline measurements

#include "chroma.h"

using namespace Chroma;

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


int main(int argc, char *argv[]) 
{
  ChromaInitialize(&argc, &argv);
  
  // Chroma Init stuff -- Open DATA and XMLDAT
  linkageHack();

  XMLWriter& xml_out = TheXMLOutputWriter::Instance();
  push(xml_out, "t_inline");

  XMLReader xml_in("./DATA");

  // Input parameter structure
  Inline_input_t  input;
  read(xml_in, "/t_inline", input);

  // Write out the input
  write(xml_out, "Input", xml_in);

  Layout::setLattSize(input.param.nrow);
  Layout::create();

  proginfo(xml_out);    // Print out basic program info

  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    
    gaugeStartup(file_xml, config_xml, u, input.cfg);

    // Write out the config header
    push(xml_out, "Config_info");
    write(xml_out, "file_xml", file_xml);
    write(xml_out, "config_xml", config_xml);
    pop(xml_out);
  }
  
  // Get the measurements
  std::istringstream Measurements_is(input.param.inline_measurement_xml);
  XMLReader MeasXML(Measurements_is);
  multi1d < Handle< AbsInlineMeasurement > > the_measurements;
  read(MeasXML, "/InlineMeasurements", the_measurements);

  QDPIO::cout << "There are " << the_measurements.size() << " measurements " << endl;

  // Measure inline observables 
  push(xml_out, "InlineObservables");

  QDPIO::cout << "Doing " << the_measurements.size() 
	      <<" measurements" << endl;
  unsigned long cur_update = 42;
  for(int m=0; m < the_measurements.size(); m++) 
  {
    AbsInlineMeasurement& the_meas = *(the_measurements[m]);
    if( cur_update % the_meas.getFrequency() == 0 ) 
    {
      // Caller writes elem rule
      push(xml_out, "elem");
      the_meas(u, cur_update, xml_out);
      pop(xml_out); 
    }
  }
  pop(xml_out); // pop("InlineObservables");

  pop(xml_out);

  ChromaFinalize();
  exit(0);
}

