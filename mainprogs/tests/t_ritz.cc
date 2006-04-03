// $Id: t_ritz.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $

#include "chroma.h"

using namespace Chroma;


// Struct for test parameters
//
typedef struct {
  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> rng_seed;
  Cfg_t        cfg;
  
  Real quark_mass;
  Real rsd_r;
  Real rsd_a;
  Real rsd_zero;
  bool projApsiP;
  int  n_renorm;
  int  n_min;
  int  n_eig;
  int max_cg;
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
    read(reader, "/params/quarkMass", params.quark_mass);

    read(reader, "/params/Cfg", params.cfg);
   
    read(reader, "/params/eig/RsdR", params.rsd_r);
    read(reader, "/params/eig/RsdA", params.rsd_a);
    read(reader, "/params/eig/RsdZero", params.rsd_zero);
    read(reader, "/params/eig/ProjApsiP",  params.projApsiP);
    read(reader, "/params/eig/MaxCG", params.max_cg);
    read(reader, "/params/eig/Nrenorm", params.n_renorm);
    read(reader, "/params/eig/Nmin", params.n_min);
    read(reader, "/params/eig/Neig", params.n_eig);

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
  write(writer, "RsdA", params.rsd_a);
  write(writer, "RsdZero", params.rsd_zero);
  write(writer, "ProjApsiP",  params.projApsiP);
  pop(writer); // Eig

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

  // Typedefs to save typing
  typedef LatticeFermion               T;
  typedef multi1d<LatticeColorMatrix>  P;
  typedef multi1d<LatticeColorMatrix>  Q;

  // Create a FermBC
  Handle< CreateFermState<T,P,Q> >  cfs(new CreateSimpleFermState<T,P,Q>(params.boundary));
  
  // The Gauge Field
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, params.cfg);

  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.
  UnprecWilsonFermAct  S_w(cfs, params.quark_mass);

  Handle< FermState<T,P,Q> > connect_state(S_w.createState(u));

  Handle< LinearOperator<T> > MM(S_w.lMdagM(connect_state));

  // Try and get lowest eigenvalue of MM
  multi1d<Real> lambda(params.n_eig);
  multi1d<Real> check_norm(params.n_eig);
  multi1d<LatticeFermion> psi(params.n_eig);
  
  int n_CG_count;

  {
    XMLBufferWriter eig_spec_xml;

    for(int i =0; i < params.n_eig; i++)
      gaussian(psi[i]);

    EigSpecRitzCG(*MM, 
		  lambda, 
		  psi, 
		  params.n_eig,
		  params.n_renorm, 
		  params.n_min, 
		  params.max_cg,
		  params.rsd_r,
		  params.rsd_a,
		  params.rsd_zero,
		  params.projApsiP,
		  n_CG_count,
		  eig_spec_xml);

    write(xml_out,"LowestEv",eig_spec_xml);
  }


  {
    QDPIO::cout << "Look for highest ev" << endl;

    Handle< LinearOperator<T> > MinusMM(new lopscl<LatticeFermion, Real>(MM, Real(-1.0)));
  
    // Look for highest ev
    for(int i =0; i < params.n_eig; i++)
      gaussian(psi[i]);
    
    QDPIO::cout << "ritz call" << endl;

    XMLBufferWriter eig_spec_xml;

    EigSpecRitzCG(*MinusMM,
		  lambda,
		  psi,
		  params.n_eig,
		  params.n_renorm,
		  params.n_min,
		  params.max_cg,
		  params.rsd_r,
		  params.rsd_a,
		  params.rsd_zero,
		  params.projApsiP,
		  n_CG_count,
		  eig_spec_xml);

    write(xml_out,"HighestEv",eig_spec_xml);
  }

  pop(xml_out);
  Chroma::finalize();
    
  exit(0);
}
