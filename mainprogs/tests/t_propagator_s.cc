// $Id: t_propagator_s.cc,v 1.9 2004-01-07 13:50:09 bjoo Exp $
/*! \file
 *  \brief Main code for propagator generation
 */

#include <iostream>
#include <cstdio>

#define MAIN

// Include everything...
#include "chroma.h"

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

      QDPIO::cout << " PROPAGATOR: Propagator for Staggered fermions" << endl;

//      read(paramtop, "numKappa", input.param.numKappa);
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

//
//  **** Pion correlator code *****
//
// These will put in a separate file one day
//



//! Function object used for constructing the time-slice set
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


//
//
//


void mesons_s(const LatticePropagator & psi,
              const multi1d<LatticeReal>& meson_phases,
              multi2d<Real> & meson_propagator,
              const multi1d<int>& t_source,
              int j_decay)
{
  LatticeReal psi_sq;
  LatticeReal meson_tmp;


// --- pasted from mesons_w ---
// Create the time-slice set
 UnorderedSet timeslice;
 timeslice.make(TimeSliceFunc(j_decay));

 // Length of lattice in j_decay direction
 int length = timeslice.numSubsets();

 // ----- end of paste -----

  multi1d<Double> hsum(length);

  int t0;
  int t;
  int t_eff;
  int n;



  t0 = t_source[j_decay];              /* Note j_decay = 0 is not permitted! */

/* Compute Psi^dag * Psi */
  psi_sq = real(trace(psi * adj(psi)));

 for(n = 0;n  < Nd; ++n )
   {
     /* Multiply with the appropriate meson phase */
          meson_tmp = meson_phases[n] * psi_sq;
     //     meson_tmp = psi_sq;

     /* Do a slice-wise sum. */
#ifdef WHY_OH_WHY
     // not sure why recode did this
     meson_tmp[0] += meson_tmp[1];
     hsum = sumMulti(meson_tmp[0], timeslice);
#endif

     hsum = sumMulti(meson_tmp, timeslice);

     //     for(t = 0;t  < ( length); ++t )
     //  {
     //  t_eff = (t - t0 + length)% length;
     //  meson_propagator[n][t_eff] += Real(hsum[t]);
     //  }


     for(t = 0;t  < ( length); ++t )
       {
         meson_propagator[n][t] = Real(hsum[t]);
       }



   }


}


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
  XMLReader xml_in("DATA");

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
  XMLFileWriter xml_out("XMLDAT");
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

  // Calcluate plaq on the transformed field
//  MesPlq(u_trans, w_plaq, s_plaq, t_plaq, link);
//  push(xml_out, "Is this gauge invariant?");
//  Write(xml_out, w_plaq);
//  Write(xml_out, s_plaq);
//  Write(xml_out, t_plaq);
//  Write(xml_out, link);
//  pop(xml_out);

//  push(xml_out, "Gauge_Field");
//  Write(xml_out, u);
//  pop(xml_out);

  xml_out.flush();

  // Phases
  // multi1d<LatticeInteger> alpha(Nd); // KS Phases
  // multi1d<LatticeInteger> beta(Nd);  // Auxiliary phases for this work (not needed here)

  // Turn on KS phases
  // New way to do this:
  //
  // for(int mu=0; mu < Nd; mu++) { 
  //   u[mu] *= StagPhases::alpha(mu);
  // }

  // Create a fermion BC. Note, the handle is on an ABSTRACT type.
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.param.boundary));

  // NOTE: THIS METHOD IS NOW OBSOLETE
//  // Apply the boundary conditions
//  setph(input.param.boundary);
//  phfctr(u);

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
  LatticeFermion tmp1, tmp2;
  const int even_source[] = { 0, 0, 0, 0 };
  multi1d<int> t_src_even(4);
  t_src_even = even_source;
			  

  const int odd_source[] = { 1, 0, 0, 0 };
  multi1d<int> t_src_odd(4);
  t_src_odd = odd_source;

  tmp1 = zero;
  srcfil(tmp1, t_src_odd ,0, 0);
  tmp2  =  zero;

  // Apply Linop
  (*D_asqtad).evenOddLinOp(tmp2, tmp1, PLUS); 

   push(xml_out, "dslash");
   Write(xml_out, tmp1);
   Write(xml_out, tmp2);
   pop(xml_out);

   Handle<const LinearOperator<LatticeFermion> > MdagM_asqtad(S_f.lMdagM(state));

  
//  srcfil(tmp1, t_src_even, 0, 0);

//  (*MdagM_asqtad)(tmp2, tmp1, PLUS);

//  push(xml_out, "MdagM");
//  Write(xml_out, tmp1);
//  Write(xml_out, tmp2);
//  pop(xml_out);

 
						       
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

  LatticeFermion q_source, psi;

  psi = zero;   // note this is ``zero'' and not 0


/*  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    QDPIO::cout << "Inversion for Color =  " << color_source << endl;
    int spin_source = 0;

    q_source = zero ;
    srcfil(q_source, t_src_even, color_source, 0);

    // Extract a fermion source
    //PropToFerm(quark_prop_source, chi, color_source, spin_source);

    //    push(xml_out, "SOURCE");
    // Write(xml_out, q_source);
   //  Write(xml_out, chi);
    // pop(xml_out);

    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    // int n_count;

    S_f.qprop(psi, state, q_source, CG_INVERTER, 
              input.param.RsdCG, input.param.MaxCG, n_count);
    
    ncg_had += n_count;
      
    push(xml_out,"Qprop");
    write(xml_out, "Mass" , input.param.Mass);
    write(xml_out, "RsdCG", input.param.RsdCG);
    Write(xml_out, n_count);
    pop(xml_out);
*/
    /*
     * Move the solution to the appropriate components
     * of quark propagator.
     */
/*    FermToProp(psi, quark_propagator, color_source, spin_source);
//    push(xml_out, "psi");
//    Write(xml_out, psi);
//    pop(xml_out);
  }
*/
/*  LatticeFermion phi = zero;
  LatticeFermion tmp3, tmp4 = zero;
  D_asqtad->evenOddLinOp(tmp1, psi, PLUS);
  phi[rb[0]] = tmp1;
  D_asqtad->evenEvenLinOp(tmp2, psi, PLUS);
  phi[rb[0]] += tmp2;
  D_asqtad->oddEvenLinOp(tmp3, psi, PLUS);
  phi[rb[1]] = tmp3;
  D_asqtad->oddOddLinOp(tmp4, psi, PLUS);
  phi[rb[1]] += tmp4;

//  push(xml_out, "tests");
//  Write(xml_out, tmp1);
//  Write(xml_out, tmp2);
//  Write(xml_out, tmp3);
//  Write(xml_out, tmp4);
//  pop(xml_out);

  push(xml_out, "PHI");
  Write(xml_out, phi);
  pop(xml_out);
*/
  
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
//   switch (input.param.prop_type) 
//   {
//   case PROP_TYPE_SZIN:
//     writeSzinQprop(quark_propagator, input.prop.prop_file, input.param.Mass);
//   break;

//  case PROP_TYPE_SCIDAC:
//    writeQprop(prop_xml, quark_propagator, input.prop.prop_file);
//    break;
 
//   default :
//     QDP_error_exit("Propagator type is unsupported.");
//   }

  //
  //  compute the pion correlator
  //
  cout << "Computing meson spectroscopy..." << endl ;

  multi1d<LatticeReal> meson_phases(Nd) ;
  int j_decay = Nd -1;
  int length = Layout::lattSize()[j_decay] ;


  MesPhas(meson_phases,j_decay) ;


  multi2d<Real>  meson_propagator(Nd,length) ;
  multi1d<Real> meson_prop(length);

    string meson_snames[4] =
      {  "M_PS"  , "M_VT"  ,
         "M_PV"  , "M_SC"
      };

    mesons_s(quark_propagator,meson_phases,meson_propagator,
             t_src_even,j_decay) ;

    push(xml_out,"Point_Point_Staggered_Hadron");
    Write(xml_out,j_decay);

    push(xml_out,"Ontology");
    string paper_reference = "Nucl.Phys.B284:299,1987, Bowler et al." ;
    Write(xml_out,paper_reference );
    string Equation = "2.18" ;
    Write(xml_out,Equation);
    pop(xml_out) ;

    push(xml_out,"Point_Point_Staggered_Meson");
    for(int stag_oper = 0; stag_oper < 4 ; stag_oper++)
      {
        for(int tt=0 ; tt < length ; ++tt)
          meson_prop[tt] = meson_propagator[stag_oper][tt] ;

        push(xml_out,meson_snames[stag_oper]);
        Write(xml_out, meson_prop);
        pop(xml_out) ;
      }
    pop(xml_out) ; // meson

    pop(xml_out) ; // hadron


  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  exit(0);
}
