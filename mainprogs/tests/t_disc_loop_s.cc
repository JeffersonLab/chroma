// $Id: t_disc_loop_s.cc,v 1.1 2004-02-09 12:41:49 mcneile Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"
//#include "meas/hadron/hadron_s.h"
/*
 *  Here we have various temporary definitions
 */
/*
enum CfgType {
  CFG_TYPE_MILC = 0,
  CFG_TYPE_NERSC,
  CFG_TYPE_SCIDAC,
  CFG_TYPE_SZIN,
  CFG_TYPE_UNKNOWN
} ;

enum PropType {
  PROP_TYPE_SCIDAC = 2,
  PROP_TYPE_SZIN,
  PROP_TYPE_UNKNOWN
} ;

enum FermType {
  FERM_TYPE_STAGGERED,
  FERM_TYPE_UNKNOWN
};
*/


using namespace QDP;


/*
 * Input 
 */


// Parameters which must be determined from the XML input
// and written to the XML output
struct Param_t
{
  FermType     FermTypeP;
  Real         Mass;      // Staggered mass
  Real         u0;        // Tadpole Factor
 
  CfgType  cfg_type;       // storage order for stored gauge configuration
  PropType prop_type;      // storage order for stored propagator

//  enum InvType  invType;            // Inverter type
  Real RsdCG;
  int MaxCG;		   // Iteration parameters
  int Nsamples; 		  // Number of stochastic sources

  Real GFAccu, OrPara;    // Gauge fixing tolerance and over-relaxation param
  int GFMax;              // Maximum gauge fixing iterations

  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> t_srce; 
};


struct Prop_t
{
  string       source_file;
  string       prop_file;
};

struct Propagator_input_t
{
  IO_version_t     io_version;
  Param_t          param;
  Cfg_t            cfg;
  Prop_t           prop;
};


//
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

//  read(inputtop, "source_file", input.source_file);
  read(inputtop, "prop_file", input.prop_file);
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
    read(inputtop, "IO_version/version", input.io_version.version);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Currently, in the supported IO versions, there is only a small difference
  // in the inputs. So, to make code simpler, extract the common bits 

  // Read the uncommon bits first
  try
  {
    XMLReader paramtop(inputtop, "param"); // push into 'param' group

    switch (input.io_version.version) 
    {
      /**************************************************************************/
    case 1 :
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << input.io_version.version << " unsupported." << endl;
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

    {
      string ferm_type_str;
      read(paramtop, "FermTypeP", ferm_type_str);
      if (ferm_type_str == "STAGGERED") {
	input.param.FermTypeP = FERM_TYPE_STAGGERED;
      } 
    }

    // GTF NOTE: I'm going to switch on FermTypeP here because I want
    // to leave open the option of treating masses differently.
    switch (input.param.FermTypeP) {
    case FERM_TYPE_STAGGERED :

      QDPIO::cout << " PROPAGATOR: Disconnected loops for Staggered fermions" << endl;

      read(paramtop, "Mass", input.param.Mass);
      read(paramtop, "u0" , input.param.u0);

#if 0
      for (int i=0; i < input.param.numKappa; ++i) {
	if (toBool(input.param.Kappa[i] < 0.0)) {
	  QDPIO::cerr << "Unreasonable value for Kappa." << endl;
	  QDPIO::cerr << "  Kappa[" << i << "] = " << input.param.Kappa[i] << endl;
	  QDP_abort(1);
	} else {
	  QDPIO::cout << " Spectroscopy Kappa: " << input.param.Kappa[i] << endl;
	}
      }
#endif

      break;

    default :
      QDP_error_exit("Fermion type not supported\n.");
    }

    {
      string cfg_type_str;
      read(paramtop, "cfg_type", cfg_type_str);
      if (cfg_type_str == "NERSC") {
	input.param.cfg_type = CFG_TYPE_NERSC;
      } else {
	QDP_error_exit("Dont know non NERSC files yet");
      }

    }

    {
      string prop_type_str;
      read(paramtop, "prop_type", prop_type_str);
      if (prop_type_str == "SZIN") {
	input.param.prop_type = PROP_TYPE_SZIN;
      } else {
	QDP_error_exit("Dont know non SZIN files yet");
      }
    }

//    read(paramtop, "invType", input.param.invType);
//    input.param.invType = CG_INVERTER;   //need to fix this
    read(paramtop, "RsdCG", input.param.RsdCG);
    read(paramtop, "MaxCG", input.param.MaxCG);
    read(paramtop, "Nsamples", input.param.Nsamples);
    read(paramtop, "GFAccu", input.param.GFAccu);
    read(paramtop, "OrPara", input.param.OrPara);
    read(paramtop, "GFMax", input.param.GFMax);

    read(paramtop, "nrow", input.param.nrow);
    read(paramtop, "boundary", input.param.boundary);
    read(paramtop, "t_srce", input.param.t_srce);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


  // Read in the gauge configuration file name
  try
  {
    read(inputtop, "Cfg", input.cfg);
    read(inputtop, "Prop", input.prop);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}

// I cant forward declare this for some reason
// Standard Time Slicery
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}

  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}

  int dir_decay;

private:
  TimeSliceFunc() {}  // hide default constructor
};

//! Propagator generation
/*! \defgroup propagator Propagator generation
 *  \ingroup main
 *
 * Main program for propagator generation. 
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DISC_DATA");

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_NERSC :
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  // 
  //  gauge invariance test
  //  


  // gauge transformed gauge fields
//  multi1d<LatticeColorMatrix> u_trans(Nd);

  // gauge transform
//  LatticeColorMatrix v ;
  
//  gaussian(v);
//  reunit(v) ; 

//  for(int dir = 0 ; dir < Nd ; ++dir)
//    {
//      u_trans[dir] = v*u[dir]*adj(shift(v,FORWARD,dir)) ;
//    }


  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_xml;

  switch (input.param.prop_type) 
  {
  case PROP_TYPE_SZIN :
//    readSzinQprop(source_xml, quark_prop_source, input.prop.source_file);
    quark_prop_source = zero;
    break;
  default :
    QDP_error_exit("Propagator type is unsupported.");
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("DISC_XMLDAT");
  push(xml_out, "propagator");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  write(xml_out, "Source_info", source_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();

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
  pop(xml_out);

  // Fix to the coulomb gauge
//  int n_gf;
  int j_decay = Nd-1;

/*  coulGauge(u, n_gf, j_decay, input.param.GFAccu, input.param.GFMax, true, input.param.OrPara);
  QDPIO::cout << "No. of gauge fixing iterations =" << n_gf << endl;

  // Calcluate plaq on the gauge fixed field
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "Is this gauge invariant?");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);
*/
  xml_out.flush();

  // Create a fermion BC. Note, the handle is on an ABSTRACT type.
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));

  //
  // Initialize fermion action
  //
  AsqtadFermAct S_f(fbc, input.param.Mass, input.param.u0);

  // Set up a state for the current u,
  // (compute fat & triple links)
  // Use S_f.createState so that S_f can pass in u0

  Handle<const ConnectState > state(S_f.createState(u));
  Handle<const EvenOddLinearOperator<LatticeFermion> > D_asqtad(S_f.linOp(state));

  // Create a fermion to apply linop to.
//  LatticeFermion tmp1, tmp2;

//  tmp1 = zero;
//  Test walfil code

//  for(int src_ind = 0; src_ind < 8; ++src_ind){
//    walfil(tmp1, 0, 3, 0, src_ind);

//  Test the volume source code
/*      
  for(int i = 0; i < 2; ++i){
    volfil(tmp1);

    tmp2  =  zero;

    // Apply Linop
    (*D_asqtad).evenOddLinOp(tmp2, tmp1, PLUS); 

    push(xml_out, "dslash");
    Write(xml_out, tmp1);
    Write(xml_out, tmp2);
    pop(xml_out);
  }
*/
  Handle<const LinearOperator<LatticeFermion> > MdagM_asqtad(S_f.lMdagM(state));

  // Machinery to do timeslice sums with
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  //
  // Loop over the source color, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a staggered propagator is a matrix in color space
  // 
  //
  LatticePropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;

  int Nsamp = input.param.Nsamples;
  int t_length = input.param.nrow[3];

  LatticeFermion q_source, psi;
  multi1d<LatticeFermion> soln_vec(Nsamp);

  LatticeComplex TrG;
//  LatticePropagator G;

  multi2d<DComplex> loop(Nsamp, t_length);
  multi1d<DComplex> sca_loop(t_length);

  sca_loop = zero;

  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0

    // Fill the volume with random noise, gaussian for now.
    // Add Z2 later.
    gaussian(q_source);

    // Compute the solution vector for the particular source
    // int n_count;

    S_f.qprop(psi, state, q_source, CG_INVERTER, 
              input.param.RsdCG, input.param.MaxCG, n_count);
    
    ncg_had += n_count;
      
    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG" , input.param.RsdCG);
    Write(xml_out, n_count);
    pop(xml_out);

    // Store the solution vectors for now
    soln_vec[i] = psi;
    
    // TrG = phase*Tr(M^-1) = Tr(conj(q_source)*psi) for scalar
    // eta' will be added soon.
    TrG = localInnerProduct(q_source, psi);

    // Do a timeslice sum

    loop[i] = sumMulti(TrG, timeslice);    
    sca_loop += loop[i];

    push(xml_out, "LOOP");
    Write(xml_out, loop[i]);
    pop(xml_out);

  } // Nsamples

//  sca_loop = sca_loop/(Real)(Nsamp);
  multi1d<Real64> re_sc(t_length), im_sc(t_length);

  // Average over stochastic samples
  for(int t = 0; t < t_length ; ++t){
    re_sc[t] = real(sca_loop[t])/Nsamp;
    im_sc[t] = imag(sca_loop[t])/Nsamp;
    sca_loop[t] = cmplx(re_sc[t], im_sc[t]);

    QDPIO::cout <<  t << sca_loop[t] <<endl;
  }
  
  // Write out timesclice sums for "offline" analysis
  // Will also write out correlators here soon!

  push(xml_out, "AV_LOOP");
  Write(xml_out, sca_loop);
  pop(xml_out);

  // Calculate the standard deviation on the average 
  // sigma_stoc = 1/(sqrt(Nsamp-1) * SD
  // write this to standard output for now
  
  multi1d<Real64> asq(t_length), s(t_length), aa(t_length), sigma(t_length);
  multi1d<Real64> imasq(t_length), ims(t_length), imaa(t_length),
                  imsigma(t_length);

  QDPIO::cout << "Here are the standard deviations on the mean" << endl;

  s = zero;

  for(int t = 0; t < t_length; ++t){
    asq[t] = pow(re_sc[t], 2);
    imasq[t] = pow(im_sc[t], 2);

    for(int i = 0; i < Nsamp; ++i){
      s[t] += pow(real(loop[i][t]), 2);
      ims[t] += pow(imag(loop[i][t]), 2);
    }
    
    aa[t] = sqrt((s[t]/Nsamp) - asq[t]);
    sigma[t] = (1/sqrt(Nsamp-1.0))*aa[t];

    imaa[t] = sqrt((ims[t]/Nsamp) - imasq[t]);
    imsigma[t] = (1/sqrt(Nsamp-1.0))*imaa[t];

    QDPIO::cout << "Real = " << sigma[t] << " Imag = " << imsigma[t] << endl;

  } 
    

  // Instantiate XML buffer to make the propagator header
  XMLBufferWriter prop_xml;
  push(prop_xml, "propagator");

  // Write out the input
  write(prop_xml, "Input", xml_in);

  // Write out the config header
  write(prop_xml, "Config_info", gauge_xml);

  // Write out the source header
  write(prop_xml, "Source_info", source_xml);

  pop(prop_xml);


  // Save the propagator
//    switch (input.param.prop_type) 
//    {
//    case PROP_TYPE_SZIN:
//     writeSzinQprop(quark_propagator, input.prop.prop_file, input.param.Mass);
//    break;

//  case PROP_TYPE_SCIDAC:
//    writeQprop(prop_xml, quark_propagator, input.prop.prop_file);
//    break;
 
//    default :
//     QDP_error_exit("Propagator type is unsupported.");
//    }


  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
