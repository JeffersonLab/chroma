// $Id: t_ritz_KS.cc,v 3.2 2007-02-22 21:11:50 bjoo Exp $

#include "chroma.h"

using namespace Chroma;

bool linkage_hack()
{
  bool foo = true;
  // All actions
  foo &= WilsonTypeFermActsEnv::registerAll();
  return foo;
}

void RitzCode5D(Handle< LinearOperatorArray<LatticeFermion> >& MM,
		const ChromaWilsonRitz_t& input,
		XMLWriter& xml_out);

void RitzCode4D(Handle< LinearOperator<LatticeFermion> >& MM,
		const ChromaWilsonRitz_t& input,
		XMLWriter& xml_out);

void RitzCode4DHw(Handle< LinearOperator<LatticeFermion> >& MM,
		  Handle< LinearOperator<LatticeFermion> >& H,
		  const ChromaWilsonRitz_t& input,
		  XMLWriter& xml_out);

  
int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  ChromaWilsonRitz_t input;
  XMLReader xml_in(Chroma::getXMLInputFileName());

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

  QDPIO::cout << "RitzEigen" << endl;

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "RitzEigen");
  proginfo(xml_out);

  // Write out the input
  write(xml_out, "Input", xml_in);
  write(xml_out, "Config_info", gauge_xml);

  xml_out.flush();

  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);
  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  // Initialise Fermion action
  std::istringstream xml_fermact_string(input.fermact);
  XMLReader fermacttop(xml_fermact_string);
  const string fermact_path = "/FermionAction";
  string fermact;
  try
  {
    read(fermacttop, fermact_path + "/FermAct", fermact);
  }
  catch (const std::string& e) 
  {
    QDPIO::cerr << "Error reading fermact: " << e << endl;
    throw;
  }

  QDPIO::cout << "FermAct = " << fermact << endl;

  // Make a reader for the stateInfo
  std::istringstream state_info_is(input.state_info);
  XMLReader state_info_xml(state_info_is);
  string state_info_path="/StateInfo";

  bool success = false;

  // Typedefs to save typing
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;

#if 1
  if( ! success ) { 
    try { 
      QDPIO::cout << "Trying 5D actions" << endl;

      // DWF-like 5D Wilson-Type stuff
      Handle< WilsonTypeFermAct5D<T,P,Q> >
	S_f(TheWilsonTypeFermAct5DFactory::Instance().createObject(fermact,
								   fermacttop,
								   fermact_path));

      Handle< FermState<T,P,Q> > state(S_f->createState(u,
							state_info_xml,
							state_info_path));

      Handle< LinearOperatorArray<LatticeFermion> > MM(S_f->lMdagM(state));


      QDPIO::cout << "Call 5D ritz code" << endl;
      RitzCode5D(MM, input, xml_out);
      QDPIO::cout << "Done with 5D ritz code" << endl;

      success = true;
    }
    catch(const std::string& e ) { 
       QDPIO::cout << "5d: " << e << endl;
    }
  }
#endif

  if( ! success ) { 
    try { 
      
      // Special case UNPRECONDITIONED_WILSON
      QDPIO::cout << "Trying 4D Wilson Like actions: " << endl;

      if( fermact == "UNPRECONDITIONED_WILSON" 
	  || fermact == "UNPRECONDITIONED_DWFTRANSF" ) {

	QDPIO::cout << "Special case. Computing Hw e-values and evecs too" << endl;
	// DWF-like 5D Wilson-Type stuff
	Handle< WilsonTypeFermAct<T,P,Q> >
	  S_f(TheWilsonTypeFermActFactory::Instance().createObject(fermact,
								   fermacttop,
								   fermact_path));
	
	Handle< FermState<T,P,Q> > state(S_f->createState(u,
							  state_info_xml,
							  state_info_path));

	Handle< LinearOperator<T> > MM(S_f->lMdagM(state));

	Handle< LinearOperator<T> > H(S_f->hermitianLinOp(state));

	RitzCode4DHw(MM, H, input, xml_out);

	success = true;
      }
      else {

	Handle< WilsonTypeFermAct<T,P,Q> >
	  S_f(TheWilsonTypeFermActFactory::Instance().createObject(fermact,
								   fermacttop,
								   fermact_path));
	
	Handle< FermState<T,P,Q> > state(S_f->createState(u,
							  state_info_xml,
							  state_info_path));

	Handle< LinearOperator<T> > MM(S_f->lMdagM(state));

	RitzCode4D(MM, input, xml_out);

	success = true;
      }
    }
    catch(const std::string& e ) { 
       QDPIO::cout << "4D: " << e << endl;
    }
  }


  pop(xml_out);
  xml_out.flush();
  xml_out.close();
  Chroma::finalize();
    
  exit(0);
}


void RitzCode4D(Handle< LinearOperator<LatticeFermion> >& MM,
		const ChromaWilsonRitz_t& input,
		XMLWriter& xml_out)
{

  // Try and get lowest eigenvalue of MM
  const Subset& s = MM->subset();
  
  multi1d<Real> lambda(input.ritz_params.Neig+input.ritz_params.Ndummy);
  multi1d<Real> check_norm(input.ritz_params.Neig);
  multi1d<LatticeFermion> psi(input.ritz_params.Neig
			      +input.ritz_params.Ndummy);
      
  
  for(int i =0; i < input.ritz_params.Neig + input.ritz_params.Ndummy; i++){
    psi[i] = zero;
    gaussian(psi[i],s);
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
  write(xml_out, "lambda_Msq", lambda); 
  
  // Check norms
  for(int i=0; i < input.ritz_params.Neig; i++) { 
    LatticeFermion Me;
    LatticeFermion lambda_e;
    (*MM)(Me, psi[i], PLUS);
    lambda_e[s] = lambda[i]*psi[i];
    
    
    LatticeFermion r_norm;
    r_norm[s] = Me - lambda_e;
    
    check_norm[i] = norm2(r_norm,s);
    check_norm[i] = sqrt(check_norm[i]);
  }
  write(xml_out, "check_norm", check_norm);
  
  for(int i=0; i < input.ritz_params.Neig; i++) {
    check_norm[i] /= fabs(lambda[i]);
  }
  write(xml_out, "check_norm_rel", check_norm);
  
  
  // Now get the absolute value of the  highest e-value
  // Work with H^{dag}H = M^{dag}M
  Real hi_RsdR = 1.0e-4;
  Real hi_RsdA = 1.0e-4;

  Handle< LinearOperator<LatticeFermion> > MinusMM = new lopscl<LatticeFermion, Real>(MM, Real(-1.0));
  
  multi1d<Real> lambda_high_aux(1);
  multi1d<LatticeFermion> lambda_high_vec(1);
  gaussian(lambda_high_vec[0],s);
  lambda_high_vec[0][s] /= sqrt(norm2(lambda_high_vec[0],s));

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
  
  lambda_high_aux[0] = fabs(lambda_high_aux[0]);
  QDPIO::cout << "|| lambda_hi || = " << lambda_high_aux[0]  << " hi_Rsd_r = " << hi_RsdR << endl;
  
  pop(high_xml);
  xml_out << high_xml;
  
  push(xml_out, "Highest");
  write(xml_out, "lambda_Msq_hi", lambda_high_aux[0]);
  pop(xml_out);
}

void RitzCode4DHw(Handle< LinearOperator<LatticeFermion> >& MM,
		  Handle< LinearOperator<LatticeFermion> >& H,
		  const ChromaWilsonRitz_t& input,
		  XMLWriter& xml_out)
{

  // Try and get lowest eigenvalue of MM
  const Subset& s = MM->subset();
  
  multi1d<Real> lambda(input.ritz_params.Neig+input.ritz_params.Ndummy);
  multi1d<Real> check_norm(input.ritz_params.Neig);
  multi1d<LatticeFermion> psi(input.ritz_params.Neig
			      +input.ritz_params.Ndummy);
      
  
  for(int i =0; i < input.ritz_params.Neig + input.ritz_params.Ndummy; i++){
    psi[i] = zero;
    gaussian(psi[i],s);
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
  write(xml_out, "lambda_Msq", lambda); 
  
  // Check norms
  for(int i=0; i < input.ritz_params.Neig; i++) { 
    LatticeFermion Me;
    LatticeFermion lambda_e;
    (*MM)(Me, psi[i], PLUS);
    lambda_e[s] = lambda[i]*psi[i];
    
    
    LatticeFermion r_norm;
    r_norm[s] = Me - lambda_e;
    
    check_norm[i] = norm2(r_norm,s);
    check_norm[i] = sqrt(check_norm[i]);
  }
  write(xml_out, "check_norm", check_norm);
  
  for(int i=0; i < input.ritz_params.Neig; i++) {
    check_norm[i] /= fabs(lambda[i]);
  }
  write(xml_out, "check_norm_rel", check_norm);
  
  // Fix to ev-s of gamma_5 wilson...
  // Try to get one:
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
  write(xml_out, "lambda_Hw", lambda);
  write(xml_out, "n_valid", n_valid);
  write(xml_out, "valid_eig", valid_eig);
  
  for(int i=0; i < input.ritz_params.Neig; i++) { 
    LatticeFermion Me;
    (*H)(Me, psi[i], PLUS);
    
    bool zeroP = toBool( fabs(lambda[i]) <  input.ritz_params.RsdZero );
    if( zeroP ) {
      check_norm[i] = norm2(Me,s);
      check_norm[i] = sqrt(check_norm[i]);
    }
    else {
      LatticeFermion lambda_e;
      LatticeFermion r_norm;
      
      lambda_e[s] = lambda[i]*psi[i];
      r_norm[s] = Me - lambda_e;
      
      
      check_norm[i] = norm2(r_norm,s);
      check_norm[i] = sqrt(check_norm[i]);
    }
   
    QDPIO::cout << "lambda_lo[" << i << "] = " << lambda[i] << "  "; 
    QDPIO::cout << "check_norm["<<i<<"] = " << check_norm[i] << endl;
  }
  write(xml_out, "check_norm", check_norm);
  
  for(int i=0; i < input.ritz_params.Neig; i++) { 
    check_norm[i] /= fabs(lambda[i]);
    QDPIO::cout << "check_norm_rel["<< i <<"] = " << check_norm[i] << endl;
  }
  QDPIO::cout << flush ;
  write(xml_out, "check_norm_rel", check_norm);
  pop(xml_out);
  
  
  // Now get the absolute value of the  highest e-value
  // Work with H^{dag}H = M^{dag}M
  Real hi_RsdR = 1.0e-4;
  Real hi_RsdA = 1.0e-4;
  
  multi1d<Real> lambda_high_aux(1);
  multi1d<LatticeFermion> lambda_high_vec(1);
  gaussian(lambda_high_vec[0],s);
  lambda_high_vec[0][s] /= sqrt(norm2(lambda_high_vec[0],s));

  int n_cg_high;
  XMLBufferWriter high_xml;
  
  Handle< LinearOperator<LatticeFermion> > MinusMM = new lopscl<LatticeFermion, Real>(MM, Real(-1.0));
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
  
  pop(high_xml);
  xml_out << high_xml;
  
  push(xml_out, "Highest");
  write(xml_out, "lambda_hi", lambda_high_aux[0]);
  pop(xml_out);

  QDPIO::cout << "Writing low eigenvalues/vectors" << endl;
  writeEigen(input, lambda, psi, lambda_high_aux[0], QDPIO_SERIAL);

}

void RitzCode5D(Handle< LinearOperatorArray<LatticeFermion> >& MM,
		const ChromaWilsonRitz_t& input,
		XMLWriter& xml_out)
{
  // Try and get lowest eigenvalue of MM
  int N5 = MM->size();
  const Subset& s = MM->subset();
  
  multi1d<Real> lambda(input.ritz_params.Neig+input.ritz_params.Ndummy);
  multi1d<Real> check_norm(input.ritz_params.Neig);
  multi2d<LatticeFermion> psi(input.ritz_params.Neig
			      +input.ritz_params.Ndummy, N5);
      
  
  for(int i =0; i < input.ritz_params.Neig + input.ritz_params.Ndummy; i++){
    for(int n=0; n < N5; n++) { 
      psi[i][n] = zero;
      gaussian(psi[i][n],s);
    }
    lambda[i] = Real(1);
  }
  
  
  int n_CG_count;
  Real delta_cycle = Real(1);
  XMLBufferWriter eig_spec_xml;
  int n_KS_count = 0;
  int n_jacob_count = 0;
  
#if 1
  QDPIO::cout << "Call EigSpecRitzKS" << endl;

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
#else

  QDPIO::cout << "Call EigSpecRitzCG" << endl;

  EigSpecRitzCG( *MM,
		 lambda,
		 psi,
		 input.ritz_params.Neig,
		 input.ritz_params.Nrenorm,
		 input.ritz_params.MinKSIter,
		 input.ritz_params.MaxCG,
		 input.ritz_params.RsdR,
		 input.ritz_params.RsdA,
		 input.ritz_params.RsdZero,
		 input.ritz_params.ProjApsiP,
		 n_CG_count,
		 eig_spec_xml);
#endif

  // Dump output
  xml_out << eig_spec_xml;
  write(xml_out, "lambda_Msq", lambda); 
  
  // Check norms
  for(int i=0; i < input.ritz_params.Neig; i++) { 
    multi1d<LatticeFermion> Me(N5);
    multi1d<LatticeFermion> lambda_e(N5);
    (*MM)(Me, psi[i], PLUS);
    for(int n =0; n < N5; n++) { 
      lambda_e[n][s] = lambda[i]*psi[i][n];
    }
    
    multi1d<LatticeFermion> r_norm(N5);
    for(int n=0; n < N5; n++) { 
      r_norm[n][s] = Me[n] - lambda_e[n];
    }
    
    check_norm[i] = norm2(r_norm[0],s);
    for(int n=1; n < N5; n++) {
      check_norm[i] += norm2(r_norm[n],s);
    }
    
    check_norm[i] = sqrt(check_norm[i]);
  }
  write(xml_out, "check_norm", check_norm);
  
  for(int i=0; i < input.ritz_params.Neig; i++) {
    check_norm[i] /= fabs(lambda[i]);
  }
  write(xml_out, "check_norm_rel", check_norm);
  
  // Now get the absolute value of the  highest e-value
  // Work with H^{dag}H = M^{dag}M
  
  
  Real hi_RsdR = 1.0e-4;
  Real hi_RsdA = 1.0e-4;
  Handle< LinearOperatorArray<LatticeFermion> > MinusMM = new lopsclArray<LatticeFermion, Real>(MM, Real(-1.0));
  
  multi1d<Real> lambda_high_aux(1);
  multi2d<LatticeFermion> lambda_high_vec(1,N5);
  for(int n=0; n < N5; n++) {
    gaussian(lambda_high_vec[0][n],s);
    lambda_high_vec[0][n][s] /= sqrt(norm2(lambda_high_vec[0],s));
  }
  
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
  
  lambda_high_aux[0] = fabs(lambda_high_aux[0]);
  QDPIO::cout << "|| lambda_hi || = " << lambda_high_aux[0]  << " hi_Rsd_r = " << hi_RsdR << endl;
  
  pop(high_xml);
  xml_out << high_xml;
  
  push(xml_out, "Highest");
  write(xml_out, "lambda_Msq_hi", lambda_high_aux[0]);
  pop(xml_out);
}
