/*
 *
 * This test program reads in: quark propagator(s);
 *                   does    : sink smearing on quark propagator(s);
 *                   returns : quark propagator(s).
 *
 * Quark propagator read in by this program will not be changed.
 *
 * sink quark smearing: arbitrary powers of laplacian;
 *                      gaussian smearing;
 *                      displacement.
 *       
 * gauge link: either "bare" link or ape-smeared link.
 *       
 */

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace QDP;




// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  CfgType         cfg_type;   // storage order for stored gauge configuration
  PropType        prop_type;  // storage order for stored propagator

  Real DummyKappa;

  Real wvf_param;
  int  WvfIntPar;
  int  LaplacePower;
  int  disp_length;
  int  disp_dir;
  WaveStateType wave_state;
  int  sink_dir;

  Real sm_fact;
  int  sm_numb;
  int  BlkMax;
  Real BlkAccu;

  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> t_srce;

  int  j_decay;
};

struct Prop_t
{
  string       prop_in_file;
  string       prop_out_file;
};

struct Propagator_input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};



void read(XMLReader& xml, const string& path, Prop_t& input);
void read(XMLReader& xml, const string& path, Propagator_input_t& input);





int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);



  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/sink_smearing", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();


  // Useful parameters that should be read from an input file
  int length = Layout::lattSize()[input.param.j_decay]; // define the temporal direction





  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;
  case CFG_TYPE_NERSC:
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }



  multi1d<LatticeColorMatrix> u_smr(Nd);
  u_smr = u;
  for(int i=0; i < input.param.sm_numb; ++i)
  {
    multi1d<LatticeColorMatrix> u_tmp(Nd);

    for(int mu = 0; mu < Nd; ++mu)
      if ( mu != input.param.j_decay )
        APE_Smear(u_smr, u_tmp[mu], mu, 0, input.param.sm_fact, 
		  input.param.BlkAccu, input.param.BlkMax, input.param.j_decay);
      else
        u_tmp[mu] = u_smr[mu];

    u_smr = u_tmp;
  }


  // Now the lattice quark propagator, just a single one for this example, together
  // with the corresponding header

  LatticePropagator quark_propagator;
  XMLReader prop_xml;
  PropHead          header;

  switch (input.param.prop_type) 
  {
  case PROP_TYPE_SZIN :
    readSzinQprop(prop_xml, quark_propagator, input.prop.prop_in_file);
    break;
    
    //  case PROP_TYPE_NERSC:
    //    readQprop("propagator_0", quark_propagator, header);
    //    break;
  default :
    QDP_error_exit("Lattice propagator type is unsupported.");
  }
  
  
  //NmlWriter nml("propagator_1");   // output in ASCII format
  //Write(nml, quark_propagator);
  //nml.close();

	


  gausSmear(u_smr, quark_propagator, input.param.wvf_param, 
	    input.param.WvfIntPar, input.param.j_decay);
	
  laplacian(u_smr, quark_propagator, input.param.j_decay, input.param.LaplacePower);
	
  displacement(u_smr,quark_propagator, input.param.disp_length, input.param.disp_dir);
	

  switch(input.param.wave_state)
  {
  case WAVE_TYPE_S_WAVE:
    break;
  case WAVE_TYPE_P_WAVE:
    {
      LatticePropagator temp;
      temp = quark_propagator;
      D_j(u_smr, temp, quark_propagator, input.param.sink_dir);
    }
    break;
  case WAVE_TYPE_D_WAVE:   
    {
      LatticePropagator temp;
      temp = quark_propagator;
      DjDk(u_smr, temp, quark_propagator, input.param.sink_dir);
    }
    break;
  default: 
    cerr<<"invaid wave_state\n";
    break;
  } 
  






  XMLFileWriter xml_out("XMLDAT_ss");  // output data file for sink_smearing
  push(xml_out, "sink_smearing");
  proginfo(xml_out);
  write(xml_out, "Input", xml_in);
  write(xml_out, "Config_info", gauge_xml);



  

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  
  push(xml_out, "Observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);


  


  // Save the propagator
  switch (input.param.prop_type) 
  {
  case PROP_TYPE_SZIN:
    writeSzinQprop(quark_propagator, input.prop.prop_out_file, input.param.DummyKappa);
    break;
    
    //  case PROP_TYPE_SCIDAC:
    //    writeQprop(prop_xml, quark_propagator, input.prop.prop_file);
    //    break;
  default :
    QDP_error_exit("Propagator type is unsupported.");
  }



  
  write(xml_out, "Input_prop_info", prop_xml);

  // Sanity check - write out the propagator (pion) correlator
  // in the Nd-1 direction
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd-1);
  
    multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator),
                                                    phases.getSet());
  
    push(xml_out, "Prop_correlator");
    Write(xml_out, prop_corr);
  }


  xml_out.flush( );
  xml_out.close( );

  exit(0);
}





void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "prop_in_file",  input.prop_in_file );
  read(inputtop, "prop_out_file", input.prop_out_file);
}






// Reader for input parameters

void read(XMLReader& xml, const string& path, Propagator_input_t& input)
{
  XMLReader inputtop(xml, path);


  // First, read the input parameter version.  Then, if this version
  // includes 'Nc' and 'Nd', verify they agree with values compiled
  // into QDP++

  // Read in the IO_version
  try
  {
    read(inputtop, "IO_version", input.io_version);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Currently, in the supported IO versions, there is only a small difference
  // in the inputs. So, to make code simpler, extract the common bits 

  try
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    switch (input.io_version.version) 
    {
    case 2:
    case 3: break;

    default:
      QDPIO::cerr << "Input parameter version " << input.io_version.version 
		  << " unsupported." << endl;
      QDP_abort(1);
    }
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }

  // Read the common bits
  try 
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    read(paramtop, "DummyKappa",   input.param.DummyKappa  );
    read(paramtop, "wvf_param",    input.param.wvf_param   );
    read(paramtop, "WvfIntPar",    input.param.WvfIntPar   );
    read(paramtop, "LaplacePower", input.param.LaplacePower);
    read(paramtop, "disp_length",  input.param.disp_length );
    read(paramtop, "disp_dir",     input.param.disp_dir    );
    read(paramtop, "wave_state",   input.param.wave_state  );
    read(paramtop, "sink_dir",     input.param.sink_dir    );
    read(paramtop, "sm_fact",      input.param.sm_fact     );
    read(paramtop, "sm_numb",      input.param.sm_numb     );
    read(paramtop, "BlkMax",       input.param.BlkMax      );
    read(paramtop, "BlkAccu",      input.param.BlkAccu     );

    // description (same order):
    // kappa value (dummy)
    // smearing width
    // number of iteration for smearing
    // number of iteration for smearing
    // power of laplacian operator
    // displacement length
    // displacement direction: x(0),y(1),z(2)
    // wave_state: "S_WAVE", "P_WAVE", or "D_WAVE"
    // sink_dir: direction of derivative at sink
    // smearing factor
    // number of smearing hits
    // Maximum number of blocking/smearing iterations
    // Blocking/smearing accuracy

    
    read(paramtop, "cfg_type",  input.param.cfg_type );
    read(paramtop, "prop_type", input.param.prop_type);



    read(paramtop, "nrow",     input.param.nrow    );
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce",   input.param.t_srce  );
    read(paramtop, "j_decay",  input.param.j_decay );


    cout << "wvf_param = " << input.param.wvf_param << endl;
    cout << "WvfIntPar = " << input.param.WvfIntPar << endl;
    cout << "LaplacePower = " << input.param.LaplacePower << endl;
    cout << "disp_length = " << input.param.disp_length << endl;
    cout << "disp_dir = " << input.param.disp_dir << endl;
    cout << "wave_state = " << input.param.wave_state << endl;


  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg",  input.cfg );
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}
