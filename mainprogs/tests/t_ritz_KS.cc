// $Id: t_ritz_KS.cc,v 1.11 2004-05-19 11:43:19 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "io/param_io.h"
#include "io/eigen_io.h"


using namespace QDP;
using namespace std;


  
int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  ChromaWilsonRitz_t input;
  XMLReader xml_in("DATA");

  try { 
    read(xml_in, "/WilsonRitzEigen", input);
  }
  catch( const string& e ) { 
    QDPIO::cerr << "Caught Exception: " << e << endl;
    QDP_error_exit("Exiting\n");
  }

  // Setup the lattice
  Layout::setLattSize(input.nrow);
  Layout::create();


  QDP::RNG::setrn(input.seed);

  QDPIO::cout << "WilsonRitzEigen" << endl;

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    QDPIO::cout << "Reading SZIN Gauge config" << endl;
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;

  case CFG_TYPE_SZINQIO:
    QDPIO::cout << "Reading SZIN QIO gauge config" << endl;
    readGauge(gauge_file_xml, gauge_xml, u, input.cfg.cfg_file, QDPIO_SERIAL);
    break;

  case CFG_TYPE_NERSC:
    QDPIO::cout << "Reading NERSC gauge config" << endl;
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  case CFG_TYPE_DISORDERED:
    QDPIO::cout << "Starting up disordered (random/hot) config" << endl;
    for(int dim=0; dim < Nd; dim++) { 
	random(u[dim]);
	reunit(u[dim]);
    }
    break;
  case CFG_TYPE_UNIT:
    QDPIO::cout << "Starting up unit gauge (free) config" << endl;
    for(int dim=0; dim < Nd; dim++) { 
	u[dim] = Real(1);
    }
    break; 
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "WilsonRitzEigen");

  write((XMLWriter &)xml_out, "InputParams", input);

  // Measure the plaquette on the gauge
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "plaquette");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link" , link);
  pop(xml_out);

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.
  // Create a FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.boundary));
  UnprecWilsonFermAct  S_w(fbc, input.Mass);

  Handle< const ConnectState > connect_state = S_w.createState(u);
  Handle< const LinearOperator<LatticeFermion> > MM = S_w.lMdagM(connect_state);

 
  // Try and get lowest eigenvalue of MM
  multi1d<Real> lambda(input.ritz_params.Neig+input.ritz_params.Ndummy);
  multi1d<Real> check_norm(input.ritz_params.Neig);
  multi1d<LatticeFermion> psi(input.ritz_params.Neig+input.ritz_params.Ndummy);

  // Initialise evecs to noise
  for(int i =0; i < input.ritz_params.Neig + input.ritz_params.Ndummy; i++) { 
    gaussian(psi[i]);
    lambda[i] = Real(1);
  }

 
  int n_CG_count;
  Real delta_cycle = Real(1);
  XMLBufferWriter eig_spec_xml;
  int n_KS_count = 0;
  int n_jacob_count = 0;
  EigSpecRitzKS(*MM, 
		lambda, 
		psi, 
		input.ritz_params.Neig,
		input.ritz_params.Ndummy,                // No of dummies
		input.ritz_params.Nrenorm, 
		input.ritz_params.MinKSIter, 
		input.ritz_params.MaxKSIter,             // Max iters / KS cycle
		input.ritz_params.MaxKS,            // Max no of KS cycles
		input.ritz_params.GammaFactor,       // Gamma factor
		input.ritz_params.MaxCG,
		input.ritz_params.RsdR,
		input.ritz_params.RsdA,  
		input.ritz_params.RsdZero,
		input.ritz_params.ProjApsiP,
		n_CG_count,
		n_KS_count,
		n_jacob_count,
		eig_spec_xml);

  // Dump output
  xml_out << eig_spec_xml;
  write(xml_out, "lambda", lambda); 

  // Check norms
  for(int i=0; i < input.ritz_params.Neig; i++) { 

    LatticeFermion Me;
    LatticeFermion lambda_e;

    (*MM)(Me, psi[i], PLUS);

    lambda_e = lambda[i]*psi[i];

    LatticeFermion r_norm = Me - lambda_e;
    check_norm[i] = sqrt(norm2(r_norm))/fabs(lambda[i]);
  }
  write(xml_out, "check_norm_rel", check_norm);

  // Fix to ev-s of gamma_5 wilson...
  // Try to get one:
  Handle< const LinearOperator<LatticeFermion> > H = S_w.gamma5HermLinOp(connect_state);  

  multi1d<bool> valid_eig(input.ritz_params.Neig);
  int n_valid;
  int n_jacob;

  fixMMev2Mev(*H, 
	      lambda, 
	      psi, 
	      input.ritz_params.Neig, 
	      input.ritz_params.RsdR,
	      input.ritz_params.RsdA, 
	      input.ritz_params.RsdZero, 
	      valid_eig, 
	      n_valid, 
	      n_jacob);

  push(xml_out, "eigFix");
  write(xml_out, "lambda", lambda);
  write(xml_out, "n_valid", n_valid);
  write(xml_out, "valid_eig", valid_eig);

  for(int i=0; i < input.ritz_params.Neig; i++) { 
    LatticeFermion Me;
    (*H)(Me, psi[i], PLUS);

    bool zeroP = toBool( fabs(lambda[i]) < input.ritz_params.RsdZero );
    if( zeroP ) {
      check_norm[i] = sqrt(norm2(Me));
    }
    else {
      LatticeFermion lambda_e;
   
      lambda_e = lambda[i]*psi[i];
      LatticeFermion r_norm = Me - lambda_e;
      check_norm[i] = sqrt(norm2(r_norm))/fabs(lambda[i]);
    }

    QDPIO::cout << "check_norm["<<i+1<<"] = " << check_norm[i] << endl;
  }
  write(xml_out, "check_norm", check_norm);
  pop(xml_out);


  multi1d<Real> szin_lambda(lambda);

  for(int i=0; i < input.ritz_params.Neig;i++ ) {
    szin_lambda[i] /= (Nd + input.Mass);
  }
  write(xml_out, "szinLamda", szin_lambda);

 
  // Now get the absolute value of the  highest e-value
  // Work with H^{dag}H = M^{dag}M
  Handle<const LinearOperator<LatticeFermion> > MinusMM = new lopscl<LatticeFermion, Real>(MM, Real(-1.0));


  Real hi_RsdR = 1.0e-4;
  Real hi_RsdA = 1.0e-4;
  
  multi1d<Real> lambda_high_aux(1);
  multi1d<LatticeFermion> lambda_high_vec(1);
  gaussian(lambda_high_vec[0]);
  lambda_high_vec[0] /= sqrt(norm2(lambda_high_vec[0]));
  int n_cg_high;
  XMLBufferWriter high_xml;

  // Initial guess -- upper bound on spectrum
  lambda_high_aux[0] = Real(8);


  push(high_xml, "LambdaHighRitz");

  // Minus MM ought to produce a negative e-value
  // since MM is herm_pos_def
  // ie minus MM is hermitian -ve definite
  EigSpecRitzCG( *MinusMM,
		 lambda_high_aux,
		 lambda_high_vec,
		 1,
		 input.ritz_params.Nrenorm,
		 input.ritz_params.MinKSIter,
		 input.ritz_params.MaxCG,
		 hi_RsdR,
		 hi_RsdA,
		 input.ritz_params.RsdZero,
		 input.ritz_params.ProjApsiP,
		 n_cg_high,
		 high_xml);

  lambda_high_aux[0] = sqrt(fabs(lambda_high_aux[0]));
  QDPIO::cout << "|| lambda_hi || = " << lambda_high_aux[0]  << " hi_Rsd_r = " << hi_RsdR << endl;

  xml_out << high_xml;

  push(xml_out, "Highest");
  write(xml_out, "lambda_hi", lambda_high_aux[0]);
  write(xml_out, "lambda_hi_szin", Real(lambda_high_aux[0]/(Real(Nd) + input.Mass)));
  pop(xml_out);

  QDPIO::cout << "Writing low eigenvalues/vectors" << endl;
  writeEigen(input, lambda, psi, lambda_high_aux[0], QDPIO_SERIAL);

  pop(xml_out);
  xml_out.close();
  QDP_finalize();
    
  exit(0);
}
