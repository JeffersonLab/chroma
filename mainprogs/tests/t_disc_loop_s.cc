// $Id: t_disc_loop_s.cc,v 1.5 2004-11-20 19:28:48 mcneile Exp $
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

  InvertParam_t  invParam;

  int Nsamples; 		  // Number of stochastic sources
  int CFGNO;		  // Configuration Number used for seeding rng
                          // This WILL be changed soon 
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
  int version;
  try
  {
    read(inputtop, "IO_version/version", version);
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

    switch (version) 
    {
      /**************************************************************************/
    case 2 :
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
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
      string prop_type_str;
      read(paramtop, "prop_type", prop_type_str);
      if (prop_type_str == "SZIN") {
	input.param.prop_type = PROP_TYPE_SZIN;
      } else {
	QDP_error_exit("Dont know non SZIN files yet");
      }
    }

//    read(paramtop, "invType", input.param.invType);
    input.param.invParam.invType = CG_INVERTER;   //need to fix this
    read(paramtop, "RsdCG", input.param.invParam.RsdCG);
    read(paramtop, "MaxCG", input.param.invParam.MaxCG);
    read(paramtop, "Nsamples", input.param.Nsamples);
    read(paramtop, "CFGNO", input.param.CFGNO);
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

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_NERSC :
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  // Read in the source along with relevant information.
  LatticeStaggeredPropagator quark_prop_source;
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
  push(xml_out, "DISCONNECTED");

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
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  // Fix to the coulomb gauge
  int n_gf;
  int j_decay = Nd-1;

  //  coulGauge(u, n_gf, j_decay, input.param.GFAccu, input.param.GFMax, true, input.param.OrPara); 
  //  QDPIO::cout << "No. of gauge fixing iterations =" << n_gf << endl;

  // Calcluate plaq on the gauge fixed field
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "Is_this_gauge_invariant");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  xml_out.flush();

  // Create a fermion BC. Note, the handle is on an ABSTRACT type.
  Handle< FermBC<LatticeStaggeredFermion> >  fbc(new SimpleFermBC<LatticeStaggeredFermion>(input.param.boundary));

  //
  // Initialize fermion action
  //
  AsqtadFermAct S_f(fbc, input.param.Mass, input.param.u0);

  // Set up a state for the current u,
  // (compute fat & triple links)
  // Use S_f.createState so that S_f can pass in u0

  Handle<const ConnectState > state(S_f.createState(u));
  Handle<const EvenOddLinearOperator<LatticeStaggeredFermion> > D_asqtad(S_f.linOp(state));
  Handle<const LinearOperator<LatticeStaggeredFermion> > MdagM_asqtad(S_f.lMdagM(state));

  // Machinery to do timeslice sums with
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  //
  // Loop over the source color, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a staggered propagator is a matrix in color space
  // 
  //
  LatticeStaggeredPropagator quark_propagator;
  XMLBufferWriter xml_buf;
  int ncg_had = 0;
  int n_count;

  int Nsamp = input.param.Nsamples;
  int t_length = input.param.nrow[3];

  LatticeStaggeredFermion q_source, psi, psi_sca1, psi_eta3, psi_eta4;
  //  multi1d<LatticeFermion> soln_vec(Nsamp);

  // Timeslice sums of fermion loops.
  // G_s is local scalar, G_eta3 is 3link eta', G_eta4 is 4link eta'
  LatticeComplex TrG_s0, TrG_s1, TrG_eta3, TrG_eta4, corr_fn_s, corr_fn_p ;

  multi2d<DComplex> loop_s0(Nsamp, t_length), loop_s1(Nsamp, t_length), 
                    loop_eta3(Nsamp, t_length), loop_eta4(Nsamp, t_length);
  multi1d<DComplex> sca0_loop(t_length), sca1_loop(t_length), eta3_loop(t_length), 
                    eta4_loop(t_length), conn_corr_s(t_length), conn_corr_p(t_length);

  sca0_loop = sca1_loop = eta3_loop = eta4_loop = zero;

  using namespace StagPhases;

  // Connected Correlator, use a point source
  // THIS IS ONLY FOR THE SCALAR AT THE MOMENT.
  // MORE THOUGHT IS NEEDED FOR THE ETA' BECAUSE OF NON-LOCALITY OF OPERATOR
  psi = zero;
  for(int color_source = 0; color_source < Nc; ++color_source) {
    int spin_source = 0;
    q_source = zero;
    srcfil(q_source, input.param.t_srce, color_source);

    S_f.qprop(psi, state, q_source, input.param.invParam, n_count);

    ncg_had += n_count;

    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG", input.param.invParam.RsdCG);
    write(xml_out, "n_count", n_count);
    pop(xml_out);

    FermToProp(psi, quark_propagator, color_source);
  }

  corr_fn_s = - alpha(1)*beta(0)*trace(adj(quark_propagator)*quark_propagator);
  conn_corr_s = sumMulti(corr_fn_s, timeslice);

  corr_fn_p = trace(adj(quark_propagator)*quark_propagator);
  conn_corr_p = sumMulti(corr_fn_p, timeslice);

  push(xml_out, "CONNECTED");
  write(xml_out, "local_scalar", conn_corr_s);
  write(xml_out, "goldstone_pion", conn_corr_p);
  pop(xml_out);

  // Seed the RNG with the cfg number for now
  Seed seed;
  seed = input.param.CFGNO;
  RNG::setrn(seed);

  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0

    // Fill the volume with random noise, gaussian for now.
    // Add Z2 later.
    RNG::savern(seed);
    QDPIO::cout << "SEED = " << seed << endl;
    //    gaussian(q_source);
    z2_src(q_source);

    // Compute the solution vector for the particular source
    // int n_count;

    S_f.qprop(psi, state, q_source, input.param.invParam,
              n_count);
    
    ncg_had += n_count;
      
    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG" , input.param.invParam.RsdCG);
    write(xml_out, "n_count", n_count);
    pop(xml_out);

    // Store the solution vectors for now; scrap this idea!
    // soln_vec[i] = psi;
    
    // TrG_s = phase*Tr(M^-1) = Tr(conj(q_source)*psi) for scalar

    TrG_s0 = localInnerProduct(q_source, psi);
    
    // 1-link scalar
    
    // The non-local operators have been commented out for now because
    // of rubbish numbers coming out of tests. More thought is needed,
    // especially for the non-local time operators (1-link scalar and
    // 4-link eta').

/*    psi_sca1 = shift(psi, FORWARD, 3);
    TrG_s1 = alpha(3)*localInnerProduct(q_source, psi_sca1);
  
    psi_eta4 = shift(shift(shift(shift(psi, FORWARD, 0), FORWARD, 1), 
              FORWARD, 2), FORWARD, 3);
 
    psi_eta3 = shift(shift(shift(psi, FORWARD, 0), FORWARD, 1), FORWARD, 2);

    TrG_eta3 = -alpha(1)*alpha(2)*localInnerProduct(q_source, psi_eta3);
    TrG_eta4 = -alpha(1)*alpha(2)*alpha(3)*localInnerProduct(q_source, psi_eta4);
*/  
    // Do timeslice sums
    loop_s0[i] = sumMulti(TrG_s0, timeslice);    
    sca0_loop += loop_s0[i];
/*    loop_s1[i] = sumMulti(TrG_s1, timeslice);
    sca1_loop += loop_s1[i];
    loop_eta3[i] = sumMulti(TrG_eta3, timeslice);
    eta3_loop += loop_eta3[i];
    loop_eta4[i] = sumMulti(TrG_eta4, timeslice);
    eta4_loop += loop_eta4[i];
*/
    push(xml_out, "LOOPS");
    write(xml_out, "scalar" , loop_s0[i]);
/*    write(xml_out, "scalar1" , loop_s1[i]);
    write(xml_out, "eta3" , loop_eta3[i]);
    write(xml_out, "eta4" , loop_eta4[i]);*/
    pop(xml_out);

  } // Nsamples

/*  multi1d<Real64> re_sc0(t_length), im_sc0(t_length), re_sc1(t_length), 
                  im_sc1(t_length), re_eta3(t_length),
                  im_eta3(t_length), re_eta4(t_length), im_eta4(t_length);
*/

  multi1d<Real64> sig_sc0(t_length), imsig_sc0(t_length);

  stoch_var(sca0_loop, loop_s0, sig_sc0, imsig_sc0, t_length, Nsamp);
  
  // Write out timesclice sums for "offline" analysis
  push(xml_out, "AV_FERMION_LOOP");
  write(xml_out, "SCALAR", sca0_loop);
  pop(xml_out);

  // Write output to seperate files for easy reading into fitting code

  ofstream out("disc.out");

  for (int t = 0; t < t_length; ++t){
    out << real(sca0_loop[t]) << " " << imag(sca0_loop[t]) << " "
        << sig_sc0[t] << " " << imsig_sc0[t] << endl;
  }

  out.close();

  ofstream out1("con.out");

  for (int t = 0; t < t_length; ++t){
    out1 << real(conn_corr_s[t]) <<endl;
  }

  out1.close();

  ofstream out2("gold.out");

  for (int t = 0; t < t_length; ++t){
    out2 << real(conn_corr_p[t]) << endl;
  }

  out2.close();

  pop(xml_out);

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
