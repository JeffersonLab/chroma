// $Id: t_ovlap5d_bj.cc,v 1.8 2004-05-19 11:43:19 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"

//#include "state.h"
//#include "actions/ferm/linop/lovlapms_w.h"
//#include "actions/ferm/fermacts/zolotarev_state.h"
//#include "actions/ferm/fermacts/zolotarev4d_fermact_bj_w.h"
//#include "actions/ferm/linop/lovlapms_w.h"
//#include "meas/eig/eig_w.h"
//#include "meas/hadron/srcfil.h"
//#include "actions/ferm/invert/invcg1.h"
//#include "util/ft/sftmom.h"

using namespace QDP;
using namespace std;

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
  Real  quark_mass;
  int  approx_order;
  double approx_min;
  double approx_max;
  Real rsd_cg;
  
  int max_cg;
  multi1d<Real> szin_pion;

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
   read(reader, "/params/qprop/RsdCG", params.rsd_cg);
   read(reader, "/params/qprop/MaxCG", params.max_cg);


   if( reader.count("/params/checking/szinPion") == 1 ) {
     read(reader, "/params/checking/szinPion", params.szin_pion);
     if( params.szin_pion.size() != params.nrow[3] ) { 
       QDP_error_exit("Comparison pion has wrong no of timeslices. Nrow[3] = %d szin_pion.size() = %d", params.nrow[3], params.szin_pion.size());
     }

   }
   else { 
     params.szin_pion.resize(0);
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
  pop(writer); // zolotarev

  push(writer, "qprop");
  write(writer, "RsdCG", params.rsd_cg);
  write(writer, "MaxCG", params.max_cg);
  pop(writer);

  if( params.szin_pion.size() > 0 ) { 
    push(writer, "checking");
    write(writer, "szinPion", params.szin_pion);
    pop(writer);
  }

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
		   XMLWriter& xml_out)
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
    filename << "eigenvector_" << setw(3) << setfill('0') << i;

    
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
    write(xml_out, "szin_enorms", szin_enorms);
    pop(xml_out);
  }
  pop(xml_out); // eigenvector test
}
  
int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);
  int G5 = Ns*Ns - 1;
  // Read the parameters 
  Param_t params;

  try { 
    readParams("./DATA", params);
  }
  catch(const string& s) { 
    QDPIO::cerr << "Caught exception " << s << endl;
    exit(1);
  }


  // Setup the lattice
  Layout::setLattSize(params.nrow);
  Layout::create();

  // Write out the params
  XMLFileWriter xml_out("t_ovlap5d.xml");
  push(xml_out, "overlapTest");

  dumpParams(xml_out, params);


  
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
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "plaquette");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.

  // Create a FermBC
  Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(params.boundary));
 
  Handle< FermBC< multi1d<LatticeFermion> > >  fbc5(new SimpleFermBC<multi1d<LatticeFermion> >(params.boundary));
  
  Handle<UnprecWilsonTypeFermAct<LatticeFermion> >  S_w(new UnprecWilsonFermAct(fbc, params.wilson_mass));
  
  XMLBufferWriter my_writer;

  //! N order Zolo approx, with wilson action.
  Zolotarev5DFermActArray   S(fbc5,
			      S_w, 
			      params.quark_mass,
			      params.approx_order, 
			      my_writer);
  
   const ConnectState* connect_state_ptr;
  multi1d<LatticeFermion> eigen_vecs;

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
		    xml_out);
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
  
  // Make me a linop (this callls the initialise function)
  // Handle<const LinearOperator< multi1d< LatticeFermion > > > D_op(S.lnonHermLinOp(connect_state));

  Handle<const LinearOperator< multi1d< LatticeFermion > > > D_op(S.lnonHermLinOp(connect_state));

  //  Handle<const LinearOperator< multi1d< LatticeFermion > > > D_MM(S.lMdagM(connect_state));


  push(xml_out, "Zolotarev5DInternal");
  xml_out << my_writer;
  pop(xml_out);

  int N5 = S.size();

  LatticeFermion chi4;
  LatticeFermion psi4;
  LatticeFermion tmp;
  int n_count;


  gaussian(chi4);
  chi4 /= sqrt(norm2(chi4));

 
  // We solve M_5d psi = gamma_5 chi
  multi1d<LatticeFermion> psi( N5 );
  multi1d<LatticeFermion> chi( N5 );

  for(int i = 0; i < N5; i++) { 
    psi[i] = zero;
    chi[i] = zero;
  }
  
  chi[N5-1] = Gamma(G5)*chi4;             // 4D source in last component
                                         // is gamma5 chi

  multi1d<LatticeFermion> tmp5_1( N5 );

  // CGNE
  //
  //  tmp5_1 = ( M^{dag} M )^{-1} chi
  //
  //  M^{dag} M is hermitian so M^{dag}M = M M^{dag}
  //
  // hence tmp5_1 = M^{-dag} M^{-1} chi
  
  // Put solution into psi  check inverse
  // (*D_op)(tmp5_1, chi, MINUS);
  InvBiCGStab( *D_op, chi, psi, params.rsd_cg, params.max_cg, n_count);

 
  // Multiply back to check inverse
  (*D_op)(tmp5_1, psi, PLUS);

  
  multi1d<LatticeFermion> tmp5_2( N5 );
  Double dnorm = Double(0);
  Double g5chi_norm = Double(0);
  for(int i = 0; i < N5; i++) { 
    tmp5_2[i] = chi[i] - tmp5_1[i];
    dnorm += norm2(tmp5_2[i]);
    g5chi_norm += norm2(chi[i]);
  }
  
  
  Double r_norm5 = sqrt(dnorm/g5chi_norm);
  QDPIO::cout << "|| chi - D psi ||/ || chi || = " << r_norm5 << endl;
  push(xml_out, "Inv5DCheck");
  write(xml_out, "r_norm", r_norm5);
  pop(xml_out);

  // Get the bottom component into psi4
  psi4 = psi[N5-1];

  // Scale to get right normalsation. 
  // Note this only works because we gamma_5-d the source initially
  psi4 *= Real(2)/(Real(1)-params.quark_mass);
  

  // This is the residuum in 4D
  LatticeFermion r;

  //! N order Zolo approx, with wilson action.
  Zolotarev4DFermAct   S4(fbc,
			  S_w, 
			  params.quark_mass,
			  params.approx_order, 
			  1.0e-9,
			  500,
			  my_writer);


  Handle<const LinearOperator< LatticeFermion > > D_op4(S4.linOp(connect_state));  

  (*D_op4)(r, psi4, PLUS);
  r -= chi4;             // Ungamma5-ed

  Double r_norm4 = sqrt(norm2(r)/norm2(chi4));
  QDPIO::cout << "|| chi4 - D psi4 || = "<< r_norm4 << endl;
  push(xml_out, "Inv4Dcheck");
  write(xml_out, "r_norm", r_norm4);
  pop(xml_out);
 

#if 0 
  // Qprop Test
  multi1d<int> coord(Nd);
  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;


  chi4 = zero;
  tmp = zero;

  // Point source
  srcfil(chi4, coord, 0, 0);
  psi4 = zero;
  
  S4.qprop(psi4, connect_state, chi4, CG_INVERTER, params.rsd_cg, params.max_cg, n_count);

  // Check original solution. psi4 contains 1/(1-m)[D^{-1} - 1 ]
  psi4 *= (Real(1)-params.quark_mass);
  psi4 += chi4; // add back contact term  psi4 = D^{-1} chi4
  (*D_op4)(r, psi4, PLUS);         // r = D psi4
  r -= chi4;                       // r = D psi4 - chi4

  r_norm4 = sqrt(norm2(r)/norm2(chi4));
  QDPIO::cout << "4D Qprop || r || = " << r_norm4 << endl;
  push(xml_out, "Consistency4D");
  write(xml_out, "r_norm", r_norm4);
  pop(xml_out);

  LatticeFermion psi4_2=zero;
  S.qprop(psi4_2, connect_state, chi4, CG_INVERTER, params.rsd_cg, params.max_cg, n_count);

  psi4_2 *= (Real(1)-params.quark_mass);  // Multiply out overall normalisation
  psi4_2 += chi4;                         // add back in contact term
  (*D_op4)(r, psi4_2, PLUS);         // r = D psi4
  r -= chi4;                       // r = D psi4 - chi4

  r_norm4 = sqrt(norm2(r)/norm2(chi4));
 
  QDPIO::cout << "5D Qprop || r || = " << r_norm4 << endl;
  push(xml_out, "Consistency5D");
  write(xml_out, "r_norm", r_norm4);
  pop(xml_out);


  // Remove contact_term
  psi4 -= chi4;
  psi4_2 -= chi4;
  //  Overlall normalisation
  psi4 *= Real(1)/( Real(1) - params.quark_mass);
  psi4_2 *= Real(1)/( Real(1) - params.quark_mass);

  r = psi4 - psi4_2;
  r_norm4 = sqrt(norm2(r)/norm2(chi4));

  QDPIO::cout << "|| 4D Qprop - 5D Qprop || = " << r_norm4 << endl;
  push(xml_out, "Qprop4D5Ddiff");
  write(xml_out, "r_norm", r_norm4);
  pop(xml_out);

#endif
  pop(xml_out);
  QDP_finalize();
    
  exit(0);
}
 
