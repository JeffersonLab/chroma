// $Id: t_ov_pbp.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"

using namespace Chroma;

enum GaugeStartType { HOT_START = 0, COLD_START = 1, FILE_START = 2 };
enum GaugeFormat { SZIN_GAUGE_FORMAT = 0, NERSC_GAUGE_FORMAT = 1 };



// Struct for test parameters
//
typedef struct {
  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> rng_seed;
  multi1d<Real> lambda;
  Real lambda_max;
  bool szin_eig;
  int gauge_start_type;
  int gauge_file_format;
  string gauge_filename;
  
  Real  wilson_mass;
  multi1d<Real>  quark_mass;  // Many quark masses
  int  approx_order;
  double approx_min;
  double approx_max;
  Real rsd_cg;
  Real rsd_cg_inner;
  
  int max_cg;
  int max_cg_inner;

  int num_noise; // Number of noise vectors
  Chirality src_chirality;

} Param_t;


// Declare routine to read the parameters
void readParams(const string& filename, Param_t& params)
{
  XMLReader reader(filename);

  try {
    // Read Params
    read(reader, "/params/lattice/nrow", params.nrow);
    read(reader, "/params/lattice/boundary", params.boundary);
    read(reader, "/params/RNG/seed", params.rng_seed);
    read(reader, "/params/zolotarev/wilsonMass", params.wilson_mass);

    // Read EigenValues: 
    //
    
    if ( reader.count("/params/SZINEValues") > 0 ) {

      // SZIN EValues
      read(reader, "/params/SZINEValues/lambda", params.lambda);
      read(reader, "/params/SZINEValues/lambdaMax", params.lambda_max);
      for(int i=0; i < params.lambda.size(); i++) { 
	params.lambda[i] *= Real(Nd) + params.wilson_mass;
      }
      params.lambda_max *= Real(Nd) + params.wilson_mass;
      params.szin_eig = true;
    }
    else if ( reader.count("/params/eValues") > 0 ) { 

      // Chroma EValues
      read(reader, "/params/eValues/lambda", params.lambda);
      read(reader, "/params/eValues/lambdaMax", params.lambda_max);
      params.szin_eig = false;
    }
    else {

      // No EValues
      params.lambda.resize(0);
      params.szin_eig = false;
    }

    read(reader, "/params/Cfg/startType", params.gauge_start_type);
    if( params.gauge_start_type == FILE_START ) { 
      read(reader, "/params/Cfg/gaugeFilename", params.gauge_filename);
      read(reader, "/params/Cfg/gaugeFileFormat", params.gauge_file_format);
    }
   

   read(reader, "/params/zolotarev/approxOrder", params.approx_order);

   if( reader.count("/params/zolotarev/approxMin") == 1 ) {
     read(reader, "/params/zolotarev/approxMin", params.approx_min);
   }
   else { 
     params.approx_min  = -1;
   }

   if( reader.count("/params/zolotarev/approxMax") == 1 ) { 
     read(reader, "/params/zolotarev/approxMax", params.approx_max);
   } else {
     params.approx_max = -1;
   }


   read(reader, "/params/zolotarev/quarkMass",  params.quark_mass);
   read(reader, "/params/zolotarev/RsdCG", params.rsd_cg_inner);
   read(reader, "/params/zolotarev/MaxCG", params.max_cg_inner);

   read(reader, "/params/PsiBarPsi/RsdCG", params.rsd_cg);
   read(reader, "/params/PsiBarPsi/MaxCG", params.max_cg);

   read(reader, "/params/PsiBarPsi/NumNoise", params.num_noise);
   int s_ch;
   read(reader,  "/params/PsiBarPsi/SrcChirality", s_ch);
 
   switch(s_ch) { 
   case 1 :
     params.src_chirality = CH_PLUS;
     break;
   case -1:
     params.src_chirality = CH_MINUS;
     break;
   case 0:
     params.src_chirality = CH_NONE;
     break;
   default:
     cerr << "Unknown value for SrcChirality " << s_ch << endl;
     exit(1);
   }


  }
  catch(const string& e) { 
    throw e;
  }
}

void dumpParams(XMLWriter& writer, Param_t& params)
{
  push(writer, "params");
  push(writer, "lattice");
  write(writer, "nrow", params.nrow);
  write(writer, "boundary", params.boundary);
  pop(writer); // lattice
  push(writer, "RNG");
  write(writer, "seed", params.rng_seed);
  pop(writer); // RNG

  if( params.lambda.size() >  0 ) { 
    if ( params.szin_eig ) {

      push(writer, "eValues");
      write(writer, "lambda", params.lambda);
      write(writer, "lambdaMax", params.lambda_max);
      pop(writer); // eValues

      push(writer, "SZINEValues");
      multi1d<Real> szin_lambda(params.lambda);
      Real szin_lambda_max = params.lambda_max;

      for(int i = 0; i < params.lambda.size(); i++) { 
	szin_lambda[i] /= (Real(Nd) + params.wilson_mass);
      }
      szin_lambda_max /= (Real(Nd) + params.wilson_mass);

      write(writer, "lambda", szin_lambda);
      write(writer, "lambdaMax", szin_lambda_max);
      pop(writer); // SZINEValues
    }
    else { 
      push(writer, "eValues");
      write(writer, "lambda", params.lambda);
      write(writer, "lambdaMax", params.lambda_max);
      pop(writer); // eValues
    }
  }
  push(writer, "Cfg");
  write(writer, "startType", params.gauge_start_type);
  if( params.gauge_start_type == FILE_START ) { 
    write(writer, "gaugeFileFormat", params.gauge_file_format);
    write(writer, "gaugeFilename", params.gauge_filename);
  }
  pop(writer); // Cfg

  push(writer, "zolotarev");
  write(writer, "approxOrder", params.approx_order);
  if( params.approx_min > 0 ) {
    write(writer, "approxMin", params.approx_min);
  }

  if( params.approx_max > 0 ) {
    write(writer, "approxMax", params.approx_max);
  }

  write(writer, "wilsonMass", params.wilson_mass);
  write(writer, "quarkMass",  params.quark_mass);
  write(writer, "RsdCG",      params.rsd_cg_inner);
  write(writer, "MaxCG",      params.max_cg_inner);
  pop(writer); // zolotarev

  push(writer, "PsiBarPsi");
  write(writer, "RsdCG", params.rsd_cg);
  write(writer, "MaxCG", params.max_cg);
  write(writer, "NumNoise", params.num_noise);

  int s_ch;
  switch( params.src_chirality ) { 
  case CH_PLUS:
    s_ch = +1;
    break;
  case CH_MINUS:
    s_ch = -1;
    break;
  case CH_NONE:
    s_ch = 0;
    break;
  default:
    cerr << "I have reached a place where I shouldnt have reached." << endl;
    exit(1);
  }
  write(writer, "SrcChirality", s_ch);

  pop(writer);

  pop(writer); // params
}

//! Read in the old SZIN eigenvectors.
//  Not only do we read the eigenvectors but we also check them
//  by computing the norm:
//
//  ||  gamma_5 D_wils e_i - lambda_i e_i ||
//
//  Since D_wils = (Nd + m_wils) D_wils_szin, we expect this
//  norm to be 
//
//   (Nd + m_wils) || gamma_5 D_wils_szin e_i - lambda_i e_i ||
//
// We also compute the old SZIN norm:
//  
//       || gamma_5 D_wils_szin e_i - lambda_i e_i || 
// 
// by dividing our original eigen nom by (Nd + m_wils)
//
// This latter norm can be checked against SZIN NMLDAT files.
//
void readEigenVecs(const multi1d<LatticeColorMatrix>& u,
		   const UnprecWilsonFermAct& S_aux,
		   const multi1d<Real>& lambda_lo, 
		   multi1d<LatticeFermion>& eigen_vec,
		   const Real wilson_mass,
		   const bool szin_eig,
		   XMLWriter& xml_out,
		   const std::string& root_prefix)
{


  // Create a connect State
  // This is where the boundary conditions are applied.

  Handle<const ConnectState>  s( S_aux.createState(u) );
  Handle<const LinearOperator<LatticeFermion> > D_w( S_aux.linOp(s) );


  // Create Space for the eigen vecs
  eigen_vec.resize(lambda_lo.size());
  
  // Create Space for the eigenvector norms
  multi1d<Real> e_norms(lambda_lo.size());
  multi1d<Real> evec_norms(lambda_lo.size());

  for(int i = 0; i < lambda_lo.size(); i++) { 
 
    // Make up the filename
    ostringstream filename;

    // this will produce eigenvector_XXX
    // where XXX is a 0 padded integer -- eg 001, 002, 010 etc
    filename << root_prefix << "eigenvector_" << setw(3) << setfill('0') << i;

    
    cout << "Reading eigenvector: " << filename.str() << endl;
    readSzinFerm(eigen_vec[i], filename.str());

    // Check e-vectors are normalized
    evec_norms[i] = (Real)sqrt(norm2(eigen_vec[i]));

    // Check the norm || Gamma_5 D e_v - lambda 
    LatticeFermion D_ev, tmp_ev, lambda_e;

    // D_ew = D ev(i)
    (*D_w)(tmp_ev, eigen_vec[i], PLUS);

    D_ev = Gamma(15)*tmp_ev;

    // Lambda_e 
    lambda_e = lambda_lo[i]*eigen_vec[i];

    D_ev -= lambda_e;

    e_norms[i] = (Real)sqrt(norm2(D_ev));        
  }    
  push(xml_out, "EigenvectorTest");
  push(xml_out, "EigenVecNorms");
  write(xml_out, "evec_norms", evec_norms);
  pop(xml_out);
  push(xml_out, "EigenTestNorms");
  write(xml_out, "e_norms", e_norms);
  pop(xml_out);

  if( szin_eig ) { 
    multi1d<Real> szin_enorms(e_norms);
    for(int i=0; i < lambda_lo.size(); i++) {
      szin_enorms[i] /= (Real(Nd)+wilson_mass);
    }
    push(xml_out, "SZINEigenTestNorms");
    write(xml_out, "szin_enorms",szin_enorms);
    pop(xml_out);
  }
  pop(xml_out); // eigenvector test
}
  
int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  std::string root_prefix="";

  if( argc == 2 ) { 
	root_prefix += argv[1];
        root_prefix += "/";
  }

  // Read the parameters 
  Param_t params;

  try { 
    readParams(root_prefix+"DATA", params);
  }
  catch(const string& s) { 
    QDPIO::cerr << "Caught exception " << s << endl;
    exit(1);
  }


  // Setup the lattice
  Layout::setLattSize(params.nrow);
  Layout::create();

  // Write out the params
  XMLFileWriter xml_out(root_prefix+"t_ov_pbp.xml");
  push(xml_out, "psiBarPsiTest");

  dumpParams(xml_out, params);


  // Create a FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(params.boundary));
  
  // The Gauge Field
  multi1d<LatticeColorMatrix> u(Nd);
  
  switch ((GaugeStartType)params.gauge_start_type) { 
  case COLD_START:
    for(int j = 0; j < Nd; j++) { 
      u(j) = Real(1);
    }
    break;
  case HOT_START:
    // Hot (disordered) start
    for(int j=0; j < Nd; j++) { 
      random(u(j));
      reunit(u(j));
    }
    break;
  case FILE_START:

    // Start from File 
    switch( (GaugeFormat)params.gauge_file_format) { 
    case SZIN_GAUGE_FORMAT:
      {
	XMLReader szin_xml;
	readSzin(szin_xml, u, params.gauge_filename);
	try { 
	  push(xml_out, "GaugeInfo");
	  xml_out << szin_xml;
	  pop(xml_out);

	}
	catch(const string& e) {
	  cerr << "Error: " << e << endl;
	}
	
      }
      break;

    case NERSC_GAUGE_FORMAT:
      {
	XMLReader nersc_xml;
	readArchiv(nersc_xml, u, params.gauge_filename);

	try { 
	  push(xml_out, "GaugeInfo");
	  xml_out << nersc_xml;
	  pop(xml_out);

	}
	catch(const string& e) {
	  cerr << "Error: " << e << endl;
	}
      }
      break;

    default:
      ostringstream file_read_error;
      file_read_error << "Unknown gauge file format" << params.gauge_file_format ;
      throw file_read_error.str();
    }
    break;
  default:
    ostringstream startup_error;
    startup_error << "Unknown start type " << params.gauge_start_type <<endl;
    throw startup_error.str();
  }


  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.
  Handle<UnprecWilsonTypeFermAct<LatticeFermion> >  S_w(new UnprecWilsonFermAct(fbc, params.wilson_mass));


  XMLBufferWriter my_writer;

  //! N order Zolo approx, with wilson action.
  // Zero quark mass explicitly, because we will use multi mass
  Zolotarev4DFermAct   S(fbc, S_w, 
			 Real(0),
			 params.approx_order, 
			 params.rsd_cg_inner,
			 params.max_cg_inner,
			 my_writer);


  const ConnectState* connect_state_ptr;
  multi1d<LatticeFermion> eigen_vecs;


  // Flick on BC's  - do not do this. Now let it be down in createState
//  phfctr(u);

  if( params.lambda.size() == 0 ) { 

    // Connect State with no eigenvectors
    connect_state_ptr = S.createState(u, 
				      Real(params.approx_min), 
				      Real(params.approx_max));
  }
  else {

    // Connect State with eigenvectors
    if( params.szin_eig ) {
      readEigenVecs(u, 
		    dynamic_cast<UnprecWilsonFermAct&>(*S_w), 
		    params.lambda, 
		    eigen_vecs, 
		    params.wilson_mass, 
		    params.szin_eig, 
		    xml_out, 
		    root_prefix);
    }
    else {
      QDP_error_exit("Non SZIN e-values not yet implmeneted");
    }

    connect_state_ptr = S.createState(u, 
				      params.lambda,
				      eigen_vecs,
				      params.lambda_max);
  }


  // Stuff the pointer into a handle. Now, the handle owns the data.
  Handle<const ConnectState> connect_state(connect_state_ptr);
						     
  // This is tricksy. I think the first size, comes into the SECOND index
  multi2d<DComplex> psi_bar_psi(params.quark_mass.size(), params.num_noise);

  for(int noise = 0; noise < params.num_noise; noise++) { 
    
    int multi_n_count=0;

    LatticeFermion noise_source;
    multi1d<LatticeFermion> inverse(params.quark_mass.size());

    // Fill source with noise
    gaussian(noise_source);

    if ( params.src_chirality != CH_NONE ) { 
      // Project with gamma_5 for positive chirality
      LatticeFermion tmp = Gamma(15)*noise_source;

      if ( params.src_chirality == CH_PLUS ) { 

	// noise source = noise_source + gamma_5 noise_source
	noise_source += tmp;
       
      }
      else { 
	noise_source -= tmp;
      }

      // noise source = (1/2) ( 1 + gamma_5 ) noise_source
      noise_source *= Real(0.5);
    }


    // Normalise noise_src to 1 (for inversion) 
    Double noise_norm = sqrt(norm2(noise_source));
    noise_source /= Real(noise_norm);


    // Fill inverse with zero (note it is "zero" not 0 )
    inverse = zero;

    
    // Invert multi mass
    int minv_ncount = 0;
  
    multi1d<Real> rsd_array(params.quark_mass.size());
    for(int m=0; m <  params.quark_mass.size(); m++) { 
      rsd_array[m] = params.rsd_cg;
    }


    // Compute 2/(1-m)[ D^{-1} - 1 ]
    S.multiQprop(inverse,
		 params.quark_mass,
		 connect_state,
		 noise_source, 
		 CG_INVERTER,
		 rsd_array,
		 params.quark_mass.size(),
		 params.max_cg,	 
		 minv_ncount);



    Real vol = Real(Layout::vol());

    for(int m=0;  m < params.quark_mass.size(); m++) {

      // Construct < source , solution > inner products for the trace
      psi_bar_psi[m][noise] = innerProduct( noise_source, inverse[m]);

      // Normalise psi-bar-psi 
      psi_bar_psi[m][noise] *= Real(2)*params.quark_mass[m];
      psi_bar_psi[m][noise] /= (Real(1) - params.quark_mass[m]*params.quark_mass[m]);
    }

  }

  for(int m=0; m < params.quark_mass.size(); m++) { 
    push(xml_out, "psiBarPsi");
    write(xml_out, "mass", params.quark_mass[m]);
    write(xml_out, "pbp", psi_bar_psi[m]);
    pop(xml_out);
  }

  pop(xml_out);
  xml_out.close();

  Chroma::finalize();
    
  exit(0);
}
