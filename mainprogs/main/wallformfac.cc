// $Id: wallformfac.cc,v 1.25 2004-06-05 21:06:27 edwards Exp $
/*! \file
 * \brief Main program for computing 3pt functions with a wall sink
 *
 * Main program for computing 3pt functions with a wall sink
 */

#include "chroma.h"

using namespace QDP;


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


struct Prop_t
{
  string       forwprop_file;
  string       backprop_file;
};


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.

  multi1d<WallFormFacType> formfac_type;
  multi1d<int> nrow;
};

struct WallFormFac_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


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




//! Main program for computing 3pt functions
/*! Main program */
int
main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  WallFormFac_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

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
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "wallFormFac");

  proginfo(xml_out);    // Print out basic program info

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config info
  write(xml_out, "Config_info", gauge_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 3);
  pop(xml_out);

  // First calculate some gauge invariant observables just for info.
  // This is really cheap.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();


  //
  // Read the quark propagator and extract headers
  //
  XMLReader forwprop_file_xml, forwprop_record_xml;
  LatticePropagator forward_quark_prop;
  ChromaProp_t forward_prop_header;
  PropSource_t forward_source_header;
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

  // Derived from input prop
  int  j_decay = forward_source_header.j_decay;
  multi1d<int> t_source = forward_source_header.t_source;

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
  PropSource_t backward_source_header;
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
   
  // Derived from input prop
  int t_sink = backward_source_header.t_source[j_decay];

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

  
  /*
   * Construct fermionic BC.
   * The BC is used to modify the U fields used for the CONSERVERED currents
   * within the formfac routines
   */
  {
    SimpleFermBC<LatticeFermion>  fbc(forward_prop_header.boundary);
    fbc.modifyU(u);   // modify the U fields
  }


  //
  // Now the 3pt contractions
  //
  SftMom phases(input.param.mom2_max, false, j_decay);

  XMLArrayWriter  xml_seq_src(xml_out, input.param.formfac_type.size());
  push(xml_seq_src, "Wilson_3Pt_fn_measurements");

  // Loop over types of form-factor
  for (int formfac_ctr = 0; formfac_ctr < input.param.formfac_type.size(); ++formfac_ctr) 
  {
    WallFormFacType formfac_type = input.param.formfac_type[formfac_ctr];
    int formfac_value = int(formfac_type);

    push(xml_seq_src);
    write(xml_seq_src, "formfac_ctr", formfac_ctr);
    write(xml_seq_src, "formfac_value", formfac_value);
    write(xml_seq_src, "formfac_type", formfac_type);

    QDPIO::cout << "Measurements for formfac_value = " << formfac_value << endl;

    WallFormFac_formfacs_t form;

    switch (formfac_type)
    {
    case WALLFF_PION:
      wallPionFormFac(form,
		      u, 
		      forward_quark_prop, backward_quark_prop, 
		      forward_quark_prop, backward_quark_prop, 
		      phases, 
		      t_source[j_decay], t_sink);
      break;

    case WALLFF_NUCL:
      wallNuclFormFac(form,
		      u, 
		      forward_quark_prop, backward_quark_prop, 
		      forward_quark_prop, backward_quark_prop, 
		      phases, 
		      t_source[j_decay], t_sink);
      break;

    case WALLFF_NUCL_CT:
    {
      /* Time-charge reverse the quark propagators */
      /* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
      LatticePropagator qf_tmp = - (Gamma(7) * forward_quark_prop * Gamma(7));
      LatticePropagator qb_tmp = - (Gamma(7) * backward_quark_prop * Gamma(7));

      wallNuclFormFac(form,
		      u, 
		      qf_tmp, qb_tmp,
		      qf_tmp, qb_tmp,
		      phases, 
		      t_source[j_decay], t_sink);
    }
    break;

    case WALLFF_DELTA:
      wallDeltaFormFac(form,
	 	       u, 
		       forward_quark_prop, backward_quark_prop, 
		       forward_quark_prop, backward_quark_prop, 
		       phases, 
		       t_source[j_decay], t_sink);
      break;

    case WALLFF_DELTA_P:
      wallDeltaPFormFac(form,
			u, 
			forward_quark_prop, backward_quark_prop, 
			forward_quark_prop, backward_quark_prop, 
			phases, 
			t_source[j_decay], t_sink);
      break;

    case WALLFF_RHO:
      wallRhoFormFac(form,
		     u, 
		     forward_quark_prop, backward_quark_prop, 
		     forward_quark_prop, backward_quark_prop, 
		     phases, 
		     t_source[j_decay], t_sink);
      break;

    case WALLFF_RHO_PI:
      wallRhoPiFormFac(form,
		       u, 
		       forward_quark_prop, backward_quark_prop, 
		       forward_quark_prop, backward_quark_prop, 
		       phases, 
		       t_source[j_decay], t_sink);
      break;

    default:
      QDPIO::cerr << "Unknown value of formfac_ctr " << formfac_ctr << endl;
      QDP_abort(1);
    }

    write(xml_seq_src, "WallFormFac", form);

    pop(xml_seq_src);   // elem
  } // end loop over formfac_ctr

  pop(xml_seq_src);  // Wilson_3Pt_fn_measurements

  // Close the output file XMLDAT
  pop(xml_out);     // wallFormFac

  xml_in.close();
  xml_out.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
