// $Id: t_ovlap_double_pass.cc,v 3.0 2006-04-03 04:59:15 edwards Exp $

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
  Real  quark_mass;
  int  approx_order;
  double approx_min;
  double approx_max;
  Real rsd_cg;
  Real rsd_cg_inner;
  
  int max_cg;
  int max_cg_inner;
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
   read(reader, "/params/zolotarev/RsdCG", params.rsd_cg_inner);
   read(reader, "/params/zolotarev/MaxCG", params.max_cg_inner);

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
  write(writer, "RsdCG",      params.rsd_cg_inner);
  write(writer, "MaxCG",      params.max_cg_inner);
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
  XMLFileWriter xml_out(root_prefix+"t_ovlap_double_pass.xml");
  push(xml_out, "overlapTest");

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
  Zolotarev4DFermAct   S(fbc, S_w, 
			   params.quark_mass,
			   params.approx_order, 
			   params.rsd_cg_inner,
			   params.max_cg_inner,
			   my_writer,
			   10,
			   OVERLAP_INNER_CG_SINGLE_PASS);


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
						     

  // Make me a linop (this callls the initialise function)
  Handle<const LinearOperator<LatticeFermion> > D_op(S.linOp(connect_state));


  // Check GW Relation.
  //
  // In massless case GWR is: D^{dag}(0) + D(0) = 2 D^{dag}(0) D(0)
  //
  // To derive for arbitrary mass case, use D(m) = m + (1-m)D(0)
  //                                    =>  D(0) = (1/m)[ D(m) - m ]
  //
  //
  // LHS is: (1/(1-m)) [ D(m) + D^{dag}(m) - 2m ] 
  //    =  (1/(1-m)) [ D(m) + D^{dag}(m) ] - 2m/(1-m)
  //
  // RHS is: (2/(1-m)^2) [ D^{dag}(m)D(m) - mD(m) - mD^{dag}(m) + m^2 ]
  //        =(1/(1-m)^2) [ 2D^{dag}(m)D(m) -2m(D(m) + D^{dag}(m)) + 2 m^2 ]
  //
  // Multiply both sides by (1-m)^2
  //
  // LHS is: (1 - m)[D(m) + D^{dag}(m) ] - 2m(1-m) 
  //       = (1 - m)[D(m) + D^{dag}(m} ] - 2m + 2m^2
  //
  // RHS is:  2D^{dag}(m)D(m) - 2m(D(m) + D^{dag}(m)) + 2 m^2
  //
  // Add 2m( D(m) + D^{dag}(m) - 2m^2 to both sides:
  //
  //  LHS is: (1 + m)[ D(m) + D^{dag}(m) ] - 2m
  //  RHS is:   2D^{dag}(m) D(m)
  //
  // 
  //  Divide both sides by 2 to get a massive version of the GWR:
  //
  //  ***********************************************************
  //
  //  ((1 + m)/2) [ D(m) + D^{dag}(m) ] - m = D^{dag}(m)D(m)
  //
  //  ***********************************************************
  //
  //  This is what we check.

  LatticeFermion psi;

  // Create some source with unit norm
  gaussian(psi);
  Double n2 = sqrt(norm2(psi));
  psi /= n2;


  LatticeFermion s1, s2, s3, tmp2;
  s1 = s2 = s3 = tmp2 = zero;

#if 0
  (*D_op)(s1,psi,PLUS);              //  s1 = D(m)  psi
  (*D_op)(s2,psi,MINUS);             //  s2 = D^{dag}(m) psi

  (*D_op)(tmp2, psi, PLUS);          //  s3 = D^{dag}(m) D(m) psi
  (*D_op)(s3, tmp2, MINUS);          //     = RHS psi

  tmp2 = s1 + s2;          // tmp2 = [ D(m) + D^{dag}(m) ] psi

 
  // tmp2 = (1 + m) [ D(m) + D^{dag}(m) ] psi
  tmp2 *= ( Real(1) + params.quark_mass )/Real(2);

  // tmp2 = (1 + m ) [ D(m) + D^{dag}(m) - m ] psi = LHS psi
  tmp2 -= params.quark_mass * psi;

  // RHS - LHS 
  s3 -= tmp2;

  // Should be zero
  Double gwr_norm = sqrt(norm2(s3));
  cout << "GWR Norm: " << gwr_norm << endl;
  write(xml_out, "gwr_norm", gwr_norm);

#endif

  //! N order Zolo approx, with wilson action.
  Zolotarev4DFermAct   S_dp(fbc, S_w, 
			    params.quark_mass,
			    params.approx_order, 
			    params.rsd_cg_inner,
			    params.max_cg_inner,
			    my_writer,
			    10,
			    OVERLAP_INNER_CG_DOUBLE_PASS);

  // Make me a linop (this callls the initialise function)
  Handle<const LinearOperator<LatticeFermion> > D_dp(S_dp.linOp(connect_state));
  
  LatticeFermion psi_dp;

  QDP::StopWatch swatch;

  swatch.reset();
  swatch.start();
  // Apply single pass
  (*D_op)(s1, psi, PLUS);
  swatch.stop();
  Double t1 = swatch.getTimeInSeconds();
  QDPIO::cout << "Single Pass algorithm took " << t1 << " seconds" << endl;
  // Apply double pass

  swatch.reset();
  swatch.start();
  (*D_dp)(s2, psi, PLUS);
  swatch.stop();
  t1 = swatch.getTimeInSeconds();
  QDPIO::cout << "Double Pass algorithm took " << t1 << " seconds" << endl;
  // Get the difference
  s3 = s2 - s1;

  Double diff = norm2(s3);

  QDPIO::cout << "Diff between single and double passes = " << sqrt(diff) << endl;


  // Now do a chiral D^{dag} D
  int G5=Ns*Ns-1;
  gaussian(psi);

  // Project onto positive chirality
  // psi <- (1/2) ( 1 + gamma_5 ) psi
  psi_dp = Gamma(G5)*psi;
  psi += psi_dp;
  psi *= 0.5;

  // Normalise
  Double norm_psi=sqrt(norm2(psi));
  psi /= Real(norm_psi);

  
  Handle< const LinearOperator<LatticeFermion> >
      MdagM_ch( S.lMdagM(connect_state, isChiralVector(psi)) );

  Handle< const LinearOperator<LatticeFermion> >
      MdagM_ch_dp( S_dp.lMdagM(connect_state, isChiralVector(psi)) );

  swatch.reset();
  swatch.start();
  (*MdagM_ch)(s1, psi, PLUS);
  swatch.stop();
  t1 = swatch.getTimeInSeconds();
  QDPIO::cout << "Single Pass algorithm took " << t1 << " seconds" << endl;


  swatch.reset();
  swatch.start();
  (*MdagM_ch_dp)(s2, psi, PLUS);
  swatch.stop();
  t1 = swatch.getTimeInSeconds();
  QDPIO::cout << "Double Pass algorithm took " << t1 << " seconds" << endl;

  s3 = s1 - s2;

  diff = norm2(s3);

  QDPIO::cout << "Diff between single and double passes = " << sqrt(diff) << endl;


  gaussian(psi);

  // Project onto positive chirality
  // psi <- (1/2) ( 1 + gamma_5 ) psi
  psi_dp = Gamma(G5)*psi;
  psi -= psi_dp;
  psi *= 0.5;

  // Normalise
  norm_psi=sqrt(norm2(psi));
  psi /= Real(norm_psi);

  Handle< const LinearOperator<LatticeFermion> >
      MdagM_ch2( S.lMdagM(connect_state, isChiralVector(psi)) );

  Handle< const LinearOperator<LatticeFermion> >
      MdagM_ch2_dp( S_dp.lMdagM(connect_state, isChiralVector(psi)) );

  swatch.reset();
  swatch.start();
  (*MdagM_ch2)(s1, psi, PLUS);
  swatch.stop();
  t1 = swatch.getTimeInSeconds();
  QDPIO::cout << "Single Pass algorithm took " << t1 << " seconds" << endl;


  swatch.reset();
  swatch.start();
  (*MdagM_ch2_dp)(s2, psi, PLUS);
  swatch.stop();
  t1 = swatch.getTimeInSeconds();
  QDPIO::cout << "Double Pass algorithm took " << t1 << " seconds" << endl;

  s3 = s1 - s2;

  diff = norm2(s3);

  QDPIO::cout << "Diff between single and double passes = " << sqrt(diff) << endl;

  
  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
