// $Id: wallformfac.cc,v 3.3 2009/11/14 20:01:46 eneil Exp $
/*! \file
 * \brief Main program for computing 3pt functions with a wall sink
 *
 * Main program for computing 3pt functions with a wall sink
 */

#include "chroma.h"

using namespace Chroma;


/*
 * Input 
 */

//! Wall-Formfactor type
enum WallFormFacType 
{
  WALLFF_PION,
  WALLFF_RHO,
  WALLFF_RHO_PI,
  WALLFF_NUCL,
  WALLFF_NUCL_CT,    // this will disappear
  WALLFF_DELTA,
  WALLFF_DELTA_P,
};


struct Prop_t
{
  string       forwprop_file;
  string       backprop_file;
};


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  int   mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
  bool  wall_source;         // use wall source or wall sink

  multi1d<WallFormFacType> formfac_type;
  multi1d<int> nrow;
};

struct WallFormFac_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


struct WallFormFac_bar_t
{
  int                     formfac_value;
  string                  formfac_type;
  WallFormFac_formfacs_t  formfacs;
};


struct WallFormFac_output_t
{
  int          out_version;

  int          mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
  bool         wall_source;         // use wall source or wall sink
  multi1d<int> nrow;

  multi1d<WallFormFac_bar_t> bar;
};



//! Read a Wall-formfactor enum
void read(XMLReader& xml, const string& path, WallFormFacType& param)
{
  string wallff_str;
  read(xml, path, wallff_str);
  if (wallff_str == "PION")
    param = WALLFF_PION;
  else if (wallff_str == "RHO")
    param = WALLFF_RHO;
  else if (wallff_str == "RHO_PI")
    param = WALLFF_RHO_PI;
  else if (wallff_str == "NUCL")
    param = WALLFF_NUCL;
  else if (wallff_str == "NUCL_CT")   // this will disappear
    param = WALLFF_NUCL_CT;
  else if (wallff_str == "DELTA")
    param = WALLFF_DELTA;
  else if (wallff_str == "DELTA_P")
    param = WALLFF_DELTA_P;
  else 
  {
    QDPIO::cerr << "Unsupported wallformfac type" << endl;
    QDP_abort(1);
  }
}


//! Write a Wall-formfactor enum
void write(XMLWriter& xml, const string& path, const WallFormFacType& param)
{
  string wallff_str;
  if (param == WALLFF_PION)
    wallff_str = "PION";
  else if (param == WALLFF_RHO)
    wallff_str = "RHO";
  else if (param == WALLFF_RHO_PI)
    wallff_str = "RHO_PI";
  else if (param == WALLFF_NUCL)
    wallff_str = "NUCL";
  else if (param == WALLFF_NUCL_CT)   // this will disappear
    wallff_str = "NUCL_CT";
  else if (param == WALLFF_DELTA)
    wallff_str = "DELTA";
  else if (param == WALLFF_DELTA_P)
    wallff_str = "DELTA_P";
  else 
  {
    QDPIO::cerr << "Unsupported formfac type" << endl;
    QDP_abort(1);
  }
  write(xml, path, wallff_str);
}


//! Propagator filenames
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "forwprop_file",input.forwprop_file);
  read(inputtop, "backprop_file",input.backprop_file);
}


//! Parameter input
void read(XMLReader& xml, const string& path, Param_t& param)
{
  XMLReader paramtop(xml, path);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 2:
    param.wall_source = false;
    break;

    /**************************************************************************/
  case 3:
    read(paramtop, "wall_source", param.wall_source);
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "Input parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "mom2_max", param.mom2_max);
  read(paramtop, "formfac_type", param.formfac_type);
  read(paramtop, "nrow", param.nrow);
}



// Reader for input parameters
void read(XMLReader& xml, const string& path, WallFormFac_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // Parameters for source construction
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the output propagator/source configuration info
    read(inputtop, "Prop", input.prop);
  }
  catch(const string& e)
  {
    QDP_error_exit("Error reading in wallformfac: %s", e.c_str());
  }
}


//! WallFormFac writer
void write(BinaryWriter& bin, const WallFormFac_bar_t& header)
{
  write(bin, header.formfac_value);
  write(bin, header.formfac_type);
  write(bin, header.formfacs);
}

//! WallFormFac writer
void write(BinaryWriter& bin, const WallFormFac_output_t& header)
{
  write(bin, header.out_version);
  write(bin, header.mom2_max);
  write(bin, header.wall_source);
  write(bin, header.nrow);
  write(bin, header.bar);
}



//! Main program for computing 3pt functions
/*! \defgroup bar3ptfnmain Computing 3pt functions
 *  \ingroup main
 *
 * Main program for computing 3pt functions
 */

int main(int argc, char *argv[])
{
  // Something breaks in this test unless Nc=3
#if QDP_NC == 3

  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();

  // Input parameter structure
  WallFormFac_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in(Chroma::getXMLInputFileName());

  // Read data
  read(xml_in, "/WallFormFac", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  QDPIO::cout << " WALLFORMFAC: Form factors for Wilson-like fermions" << endl;
  QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;
  QDPIO::cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;

  // Read in the configuration along with relevant information.
  QDPIO::cout << "Attempt to initialize the gauge field" << endl;

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  // Startup gauge
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  // Next check the gauge field configuration by reunitarizing.
  unitarityCheck(u);

  QDPIO::cout << "Gauge field successfully initialized" << endl;


  // Instantiate XML writer for XMLDAT
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

  push(xml_out, "wallFormFac");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 4);
  pop(xml_out);

  // First calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  /*
   * Read the quark propagator and extract headers
   */
  // Read the forward prop
  XMLReader forwprop_file_xml, forwprop_record_xml;
  LatticePropagator forward_quark_prop;
  ChromaProp_t forward_prop_header;
  PropSourceConst_t forward_source_header;
  {
    QDPIO::cout << "Attempt to read forward propagator" << endl;
    readQprop(forwprop_file_xml, 
	      forwprop_record_xml, forward_quark_prop,
	      input.prop.forwprop_file, QDPIO_SERIAL);
    QDPIO::cout << "Forward propagator successfully read" << endl;
   
    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try
    {
      read(forwprop_record_xml, "/Propagator/ForwardProp", forward_prop_header);
      read(forwprop_record_xml, "/Propagator/PropSource", forward_source_header);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
      throw;
    }
  }

  // Sanity check
  if (input.param.wall_source && forward_source_header.source.id != "WALL_SOURCE")
  {
    QDPIO::cerr << "Wallformfac: wall_source flag set but not a wall source forward prop" << endl;
    QDP_abort(1);
  }

  // Derived from input prop
  int  j_decay  = forward_source_header.j_decay;
  int  t_source = forward_source_header.t_source;

  // Sanity check - write out the norm2 of the forward prop in the j_decay direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    multi1d<Double> forward_prop_corr = sumMulti(localNorm2(forward_quark_prop), 
						 phases.getSet());

    push(xml_out, "Forward_prop_correlator");
    write(xml_out, "forward_prop_corr", forward_prop_corr);
    pop(xml_out);
  }

  // Save forward prop input
  push(xml_out, "ForwardPropHeaders");
  write(xml_out, "ForwardProp", forward_prop_header);
  write(xml_out, "PropSource", forward_source_header);
  pop(xml_out);


  // Read the backward propagator
  XMLReader backprop_file_xml, backprop_record_xml;
  LatticePropagator backward_quark_prop;
  ChromaProp_t backward_prop_header;
  PropSourceConst_t backward_source_header;
  {
    QDPIO::cout << "Attempt to read backward propagator" << endl;
    readQprop(backprop_file_xml, 
	      backprop_record_xml, backward_quark_prop,
	      input.prop.backprop_file, QDPIO_SERIAL);
   
    // Try to invert this record XML into a ChromaProp struct
    // Also pull out the id of this source
    try
    {
      read(backprop_record_xml, "/Propagator/ForwardProp", backward_prop_header);
      read(backprop_record_xml, "/Propagator/PropSource", backward_source_header);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << "Error extracting backward_prop header: " << e << endl;
      throw;
    }
  }
  QDPIO::cout << "Backward propagator successfully read" << endl;
   
  // Sanity check
  if (! input.param.wall_source && backward_source_header.source.id != "WALL_SOURCE")
  {
    QDPIO::cerr << "Wallformfac: wall_source flag false but not a wall source backward prop" << endl;
    QDP_abort(1);
  }

  // Derived from input prop
  int t_sink = backward_source_header.t_source;

  // Sanity check - write out the norm2 of the backward prop in the j_decay direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, j_decay);

    multi1d<Double> backward_prop_corr = sumMulti(localNorm2(backward_quark_prop), 
						 phases.getSet());

    push(xml_out, "Backward_prop_correlator");
    write(xml_out, "backward_prop_corr", backward_prop_corr);
    pop(xml_out);
  }

  // Save backward prop input
  push(xml_out, "BackwardPropHeaders");
  write(xml_out, "ForwardProp", backward_prop_header);
  write(xml_out, "PropSource", backward_source_header);
  pop(xml_out);

  
#if 0
  /*
   * Construct fermionic BC.
   * The BC is used to modify the U fields used for the CONSERVERED currents
   * within the formfac routines
   */
  {
    SimpleFermBC<LatticeFermion>  fbc(forward_prop_header.boundary);
    fbc.modifyU(u);   // modify the U fields
  }
#endif


  // Phase factors
  SftMom phases(input.param.mom2_max, false, j_decay);


  /*
   * Construct the forward propagator evaluated at the sink.
   * In the case of wall-sinks, simply wall-sink smear. This is a slice-sum.
   * For a wall-source, then use the smearing of the backward 
   * sink propagator and apply it to the forward propagator. In this
   * latter case, one evaluates the smeared prop at the sink location.
   */
  Propagator forward_quark_x2;

  if (input.param.wall_source)
  {
    // Source is a wall, so sink smear according to what was
    // used for the backward prop
    LatticePropagator forward_quark_tmp = forward_quark_prop;

    // Sink smear the propagator
    try
    {
      std::istringstream  xml_s(backward_source_header.source.xml);
      XMLReader  sinktop(xml_s);
      QDPIO::cout << "Source = " << backward_source_header.source.id << endl;

      Handle< QuarkSourceSink<LatticePropagator> >
	sinkSmearing(ThePropSinkSmearingFactory::Instance().createObject(backward_source_header.source.id,
									 sinktop,
									 backward_source_header.source.path,
									 u));
      (*sinkSmearing)(forward_quark_tmp);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "wallformfac: Caught Exception creating sink smear: " << e << endl;
      QDP_abort(1);
    }

#if 0
    // OLD code - just for reference
    if (backward_source_header.source.id == "SHELL_SOURCE")
    {
      sink_smear2(u, forward_quark_tmp, 
		  backward_source_header.sourceSmearParam.wvf_kind, 
		  backward_source_header.sourceSmearParam.wvf_param, 
		  backward_source_header.sourceSmearParam.wvfIntPar, 
		  j_decay);
    }
#endif

    // Grab the smeared forward prop at the sink location
    multi1d<int> t_snk = backward_source_header.getTSrce();
    forward_quark_x2 = peekSite(forward_quark_tmp, t_snk);
  }
  else
  {
    // Sink is a wall
    // Project forward propagator onto zero momentum: Do a slice-wise sum.
    forward_quark_x2 = sum(forward_quark_prop, phases.getSet()[t_sink]);
  }
    

  //
  // Big nested structure that is image of entire file
  //
  WallFormFac_output_t  form;
  form.bar.resize(input.param.formfac_type.size());

  form.out_version = 4;  // bump this up everytime something changes
  form.nrow = input.param.nrow;
  form.mom2_max = input.param.mom2_max;
  form.wall_source = input.param.wall_source;

  multi1d<string> wallformfac_names(7);
  wallformfac_names[0] = "PION";
  wallformfac_names[1] = "RHO";
  wallformfac_names[2] = "RHO_PI";
  wallformfac_names[3] = "NUCL";
  wallformfac_names[4] = "NUCL_CT";
  wallformfac_names[5] = "DELTA";
  wallformfac_names[6] = "DELTA_P";


  //
  // Now the 3pt contractions
  //
  XMLArrayWriter  xml_seq_src(xml_out, input.param.formfac_type.size());
  push(xml_seq_src, "Wilson_3Pt_fn_measurements");

  QDPIO::cout << "Looping over " << input.param.formfac_type.size() 
	      << " kinds of form-factors" << endl;

  // Loop over types of form-factor
  for (int formfac_ctr = 0; formfac_ctr < input.param.formfac_type.size(); ++formfac_ctr) 
  {
    WallFormFacType formfac_type = input.param.formfac_type[formfac_ctr];
    int formfac_value = int(formfac_type);

    push(xml_seq_src);
    write(xml_seq_src, "formfac_ctr", formfac_ctr);
    write(xml_seq_src, "formfac_value", formfac_value);
    write(xml_seq_src, "formfac_type", formfac_type);

    form.bar[formfac_ctr].formfac_value = formfac_value;
    form.bar[formfac_ctr].formfac_type = wallformfac_names[formfac_value];

    QDPIO::cout << "Measurements for formfac_value = " << formfac_value << endl;

    switch (formfac_type)
    {
    case WALLFF_PION:
      wallPionFormFac(form.bar[formfac_ctr].formfacs,
		      u, 
		      forward_quark_prop, backward_quark_prop, 
		      forward_quark_prop, backward_quark_prop, 
		      forward_quark_x2, forward_quark_x2,
		      phases, 
		      t_source, input.param.wall_source);
      break;

    case WALLFF_NUCL:
      wallNuclFormFac(form.bar[formfac_ctr].formfacs,
		      u, 
		      forward_quark_prop, backward_quark_prop, 
		      forward_quark_prop, backward_quark_prop, 
		      forward_quark_x2, forward_quark_x2,
		      phases, 
		      t_source, input.param.wall_source);
      break;

#if 0
    case WALLFF_NUCL_CT:
    {
      /* Time-charge reverse the quark propagators */
      /* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
      LatticePropagator qf_tmp = - (Gamma(7) * forward_quark_prop * Gamma(7));
      LatticePropagator qb_tmp = - (Gamma(7) * backward_quark_prop * Gamma(7));

      wallNuclFormFac(form.bar[formfac_ctr].formfacs,
		      u, 
		      qf_tmp, qb_tmp,
		      qf_tmp, qb_tmp,
		      forward_quark_x2, forward_quark_x2,
		      phases, 
		      t_source, input.param.wall_source);
    }
    break;
#endif

    case WALLFF_DELTA:
      wallDeltaFormFac(form.bar[formfac_ctr].formfacs,
	 	       u, 
		       forward_quark_prop, backward_quark_prop, 
		       forward_quark_prop, backward_quark_prop, 
		       forward_quark_x2, forward_quark_x2,
		       phases, 
		       t_source, input.param.wall_source);
      break;

    case WALLFF_DELTA_P:
      wallDeltaPFormFac(form.bar[formfac_ctr].formfacs,
			u, 
			forward_quark_prop, backward_quark_prop, 
			forward_quark_prop, backward_quark_prop, 
			forward_quark_x2, forward_quark_x2,
			phases, 
			t_source, input.param.wall_source);
      break;

    case WALLFF_RHO:
      wallRhoFormFac(form.bar[formfac_ctr].formfacs,
		     u, 
		     forward_quark_prop, backward_quark_prop, 
		     forward_quark_prop, backward_quark_prop, 
		     forward_quark_x2, forward_quark_x2,
		     phases, 
		     t_source, input.param.wall_source);
      break;

    case WALLFF_RHO_PI:
      wallRhoPiFormFac(form.bar[formfac_ctr].formfacs,
		       u, 
		       forward_quark_prop, backward_quark_prop, 
		       forward_quark_prop, backward_quark_prop, 
		       forward_quark_x2, forward_quark_x2,
		       phases, 
		       t_source, input.param.wall_source);
      break;

    default:
      QDPIO::cerr << "Unknown value of formfac_ctr " << formfac_ctr << endl;
      QDP_abort(1);
    }

    pop(xml_seq_src);   // elem
  } // end loop over formfac_ctr

  pop(xml_seq_src);  // Wilson_3Pt_fn_measurements

  // Close the output file XMLDAT
  pop(xml_out);     // wallFormFac

  // Dump binary output
  BinaryFileWriter  bin_out("wallformfac.dat");
  write(bin_out, form);
  bin_out.close();

  END_CODE();

  // Time to bolt
  Chroma::finalize();

#endif

  exit(0);
}
