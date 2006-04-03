// $Id: t_ritz5d_KS.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $

#include "chroma.h"

using namespace Chroma;

// Struct for test parameters
//
struct Param_t {
  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> rng_seed;
  Cfg_t        cfg;
  
  Real quark_mass;
  Real rsd_r;
  Real rsd_a;
  Real rsd_zero;
  Real gamma_factor;
  bool projApsiP;  
  int  n_eig;
  int n_dummy;
  int max_cg;
  int max_KS;

  Real wilson_mass;
  int approx_order;
  Real approx_min;
  Real approx_max;
  
};

// Declare routine to read the parameters
void readParams(const string& filename, Param_t& params)
{
  XMLReader reader(filename);

  try 
  {
    // Read Params
    read(reader, "/params/lattice/nrow", params.nrow);
    read(reader, "/params/lattice/boundary", params.boundary);
    read(reader, "/params/RNG/seed", params.rng_seed);
    read(reader, "/params/quarkMass", params.quark_mass);

    read(reader, "/params/Cfg", params.cfg);
   
    read(reader, "/params/eig/RsdR", params.rsd_r);
    read(reader, "/params/eig/RsdA", params.rsd_a);
    read(reader, "/params/eig/RsdZero", params.rsd_zero);
    read(reader, "/params/eig/gammaFactor", params.gamma_factor);
    read(reader, "/params/eig/ProjApsiP",  params.projApsiP);
    read(reader, "/params/eig/MaxKS", params.max_KS);
    read(reader, "/params/eig/MaxCG", params.max_cg);
    read(reader, "/params/eig/Neig", params.n_eig);
    read(reader, "/params/eig/Ndummy", params.n_dummy);
    read(reader, "/params/zolotarev/approxOrder", params.approx_order);
    read(reader, "/params/zolotarev/approxMin", params.approx_min);
    read(reader, "/params/zolotarev/approxMax", params.approx_max);
    read(reader, "/params/zolotarev/wilsonMass", params.wilson_mass);

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

  write(writer, "quarkMass", params.quark_mass);
  write(writer, "Cfg", params.cfg);

  push(writer, "eig");
  write(writer, "RsdR", params.rsd_r);
  write(writer, "MaxCG", params.max_cg);
  write(writer, "Neig", params.n_eig);
  write(writer, "Ndummy", params.n_dummy);
  write(writer, "RsdA", params.rsd_a);
  write(writer, "RsdZero", params.rsd_zero);
  write(writer, "ProjApsiP",  params.projApsiP);
  write(writer, "gammaParam", params.gamma_factor);
  pop(writer); // Eig

  push(writer, "zolotarev");
  write(writer, "approxOrder", params.approx_order);
  write(writer, "approxMin", params.approx_min);
  write(writer, "approxMax", params.approx_max);
  write(writer, "wilsonMass", params.wilson_mass);
  pop(writer); // Zolotarev

  pop(writer); // params
}

  
int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Read the parameters 
  Param_t params;

  try { 
    readParams(Chroma::getXMLInputFileName(), params);
  }
  catch(const string& s) { 
    QDPIO::cerr << "Caught exception " << s << endl;
    exit(1);
  }


  // Setup the lattice
  Layout::setLattSize(params.nrow);
  Layout::create();

  // Write out the params
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "ritzTest");

  dumpParams(xml_out, params);

  // The Gauge Field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, params.cfg);

  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  //! Wilsoniums;

  // Create a FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(params.boundary));

  // Auxiliary action
  Handle<UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > S_w(new UnprecWilsonFermAct(fbc, params.wilson_mass));

  Handle< FermBC< multi1d<LatticeFermion> > >  fbc5(new SimpleFermBC<multi1d<LatticeFermion> >(params.boundary));

  XMLBufferWriter my_writer;
UnprecOvlapContFrac5DFermActArray S(fbc5,
				    S_w, 
				    params.quark_mass,
				    params.approx_order,
				    my_writer);
  
  Handle< const ConnectState > connect_state(S.createState(u, 
							   Real(params.approx_min), 
							   Real(params.approx_max)));

  Handle<const LinearOperator< multi1d< LatticeFermion > > > MM(S.lMdagM(connect_state));

  // Dump Zolo Info
  xml_out << my_writer;

  // Get back 5th dim length
  int N5 = S.size();

  int n_dummy = 2;
  // Try and get lowest eigenvalue of MM
  multi1d<Real> lambda(params.n_eig+params.n_dummy);
  multi1d<Real> check_norm(params.n_eig);
  multi2d<LatticeFermion> psi(params.n_eig+params.n_dummy, N5);
  
  int n_renorm = 10;
  int n_min = 5;
  bool ProjApsiP = true;
  int n_CG_count;


  Real delta_cycle = Real(1);
  Real gamma_factor = Real(1);


  XMLBufferWriter eig_spec_xml;

  for(int i =0; i < params.n_eig+ params.n_dummy; i++) { 
    for(int n=0; n < N5; n++) { 
      gaussian(psi[i][n]);
    }
    lambda[i] = Real(1);
  }

  int n_KS_count = 0;
  int n_jacob_count = 0;
  EigSpecRitzKS(*MM, 
		lambda, 
		psi, 
		params.n_eig,
		params.n_dummy,                // No of dummies
		n_renorm, 
		n_min, 
		200,             // Max iters / KS cycle
		params.max_KS,            // Max no of KS cycles
		params.gamma_factor,       // Gamma factor
		params.max_cg,
		params.rsd_r,
		params.rsd_a,  
		params.rsd_zero,
		params.projApsiP,
		n_CG_count,
		n_KS_count,
		n_jacob_count,
		eig_spec_xml);

  xml_out << eig_spec_xml;
  write(xml_out, "lambda", lambda); 

  // Check norms
  for(int i=0; i < params.n_eig; i++) { 
    multi1d<LatticeFermion> Me(N5);
    multi1d<LatticeFermion> lambda_e(N5);
    (*MM)(Me, psi[i], PLUS);
    for(int n =0; n < N5; n++) { 
      lambda_e[n] = lambda[i]*psi[i][n];
    }

    multi1d<LatticeFermion> r_norm(N5);
    for(int n=0; n < N5; n++) { 
      r_norm[n] = Me[n] - lambda_e[n];
    }

    check_norm[i] = norm2(r_norm[0]);
    for(int n=1; n < N5; n++) {
      check_norm[i] += norm2(r_norm[n]);
    }

    check_norm[i] = sqrt(check_norm[i]);
  }
  write(xml_out, "check_norm", check_norm);

  for(int i=0; i < params.n_eig; i++) {
    check_norm[i] /= fabs(lambda[i]);
  }
  write(xml_out, "check_norm_rel", check_norm);


  // Fix to ev-s of Matrix
  // Try to get one:
  Handle< const LinearOperator< multi1d<LatticeFermion> > > M = S.linOp(connect_state);  

  multi1d<bool> valid_eig(params.n_eig);
  int n_valid;
  int n_jacob;

  fixMMev2Mev(*M, lambda, psi, params.n_eig, params.rsd_r,
	      params.rsd_a, params.rsd_zero, valid_eig, n_valid, n_jacob);

	       

  push(xml_out, "eigFix");
  write(xml_out, "lambda", lambda);
  write(xml_out, "n_valid", n_valid);
  write(xml_out, "valid_eig", valid_eig);
  for(int i=0; i < params.n_eig; i++) { 
    multi1d<LatticeFermion> Me(N5);
    (*M)(Me, psi[i], PLUS);

    bool zeroP = toBool( fabs(lambda[i]) < params.rsd_zero );
    if( zeroP ) {
      check_norm[i] = norm2(Me[0]);
      for(int n=1; n < N5; n++) { 
	check_norm[i] += norm2(Me[n]);
      }
      check_norm[i] = sqrt(check_norm[i]);
    }
    else {
      multi1d<LatticeFermion> lambda_e(N5);
      multi1d<LatticeFermion> r_norm(N5);

      for( int n=0; n < N5; n++) { 
	lambda_e[n] = lambda[i]*psi[i][n];
	r_norm[n] = Me[n] - lambda_e[n];
      }

      check_norm[i] = norm2(r_norm[0]);
      for( int n=1; n < N5; n++) { 
	check_norm[i] += norm2(r_norm[n]);
      }

      check_norm[i] = sqrt(check_norm[i]);
    }

    QDPIO::cout << "check_norm["<<i+1<<"] = " << check_norm[i] << endl;
  }
  write(xml_out, "check_norm", check_norm);

  for(int i=0; i < params.n_eig; i++) { 
    check_norm[i] /= fabs(lambda[i]);
    QDPIO::cout << "check_norm_rel["<< i+1 <<"] = " << check_norm[i] << endl;
  }
  write(xml_out, "check_norm_rel", check_norm);
  pop(xml_out);


  for(int i=0; i < params.n_eig;i++ ) {
    lambda[i] /= (Nd + params.quark_mass);
  }
  write(xml_out, "szinLamda", lambda);
  
  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
