#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "handle.h"
#include "eoprec_linop.h"
#include "eoprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"


using namespace Chroma;


struct InputParam {
  multi1d<int> nrow;
  GroupXML_t cfg;
  GroupXML_t fermact;
  GroupXML_t inv_param;
  multi1d<Real> shifts;
  int iters;
  QDP::Seed rng_seed;
};

void read(XMLReader& xml, const std::string& path, InputParam& p)
{
  XMLReader paramtop(xml, path);
  read(paramtop, "nrow", p.nrow);
  p.cfg = readXMLGroup(paramtop, "Cfg", "cfg_type");
  p.fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
  p.inv_param = readXMLGroup(paramtop, "InvertParam", "invType");
  read(paramtop, "shifts", p.shifts);
  if( paramtop.count("iters") == 0 ) { 
    p.iters = 5;
  }
  else { 
    read(paramtop, "iters", p.iters);
  }
  if (paramtop.count("RNG") > 0)
    read(paramtop, "RNG", p.rng_seed);
  else
    p.rng_seed = 11;     // default seed
}



bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= InlineAggregateEnv::registerAll();
  foo &= GaugeInitEnv::registerAll();

  return foo;
}


int main(int argc, char *argv[]) 
{


  using T = LatticeFermion;
  using Q = multi1d<LatticeColorMatrix>;
  using P = multi1d<LatticeColorMatrix>;
  using S_asymm_T = EvenOddPrecWilsonTypeFermAct<T,P,Q>;
  using LinOpAsymm_T = EvenOddPrecLinearOperator<T,P,Q>;

  Q u;
  Handle<S_asymm_T> S_asymm;
  Handle<FermState<T,P,Q> > state; 
  InputParam p;
  Handle<LinOpAsymm_T> M_asymm;

  Chroma::initialize(&argc, &argv);
  QDPIO::cout << "Linkage = " << linkageHack() << std::endl;

  // Read Params
  XMLReader xml_in;
  xml_in.open(Chroma::getXMLInputFileName());
  read(xml_in, "/ShiftedSolverTest", p);
  xml_in.close();

  // Setup the layout and stuff
  Layout::setLattSize(p.nrow);
  Layout::create();

  // Initialise the RNG
  QDP::RNG::setrn(p.rng_seed);


  // Print back some info
  u.resize(Nd); 
  StopWatch swatch;
  swatch.reset();
  swatch.start();
  try {
    XMLReader gauge_file_xml, gauge_xml;
    QDPIO::cout << "Attempt to read gauge field" << std::endl;
    std::istringstream  xml_c(p.cfg.xml);
    XMLReader  cfgtop(xml_c);
    QDPIO::cout << "Gauge initialization: cfg_type = " << p.cfg.id << std::endl;

    Handle< GaugeInit > gaugeInit(TheGaugeInitFactory::Instance().createObject(p.cfg.id, cfgtop, p.cfg.path));
    (*gaugeInit)(gauge_file_xml, gauge_xml, u);
  }
  catch(...)
  {
    // This might happen on any node, so report it
    QDPIO::cerr << "Caught generic exception during gaugeInit" << std::endl;
    throw;
  }
  swatch.stop();
  QDPIO::cout << "Gauge Field Read took " << swatch.getTimeInSeconds() << " sec" << std::endl;

  QDPIO::cout << "Reunitarizing..." ;
  for(auto dim=0; dim < Nd; ++dim) {
    QDPIO::cout << ".";
    reunit(u[dim]);
  }
  QDPIO::cout << " done" << std::endl;

  // Check plaquette for sanity of gauge read 
  {
    Double w_plaq, s_plaq, t_plaq, link;
    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    QDPIO::cout << "Plaquette = " << w_plaq <<  "( spatial =  " << s_plaq << " , temporal = " << t_plaq << " )   Link Trace = " << link << std::endl;
  }

  QDPIO::cout << "Shifted Inverter Test" << std::endl;
  QDPIO::cout << "Shifts: { " << std::endl;
  for(unsigned int j=0; j < p.shifts.size(); ++j) {
    QDPIO::cout << "\t " << p.shifts[j] << std::endl;
  }
  QDPIO::cout << "}" << std::endl;

  QDPIO::cout << "Inverter XML: " << p.inv_param.xml << std::endl;
  QDPIO::cout << "Creating Fermion Action" << std::endl;
  {
    std::istringstream fermact_xml( p.fermact.xml );
    XMLReader fermacttop(fermact_xml);

    S_asymm = dynamic_cast<S_asymm_T*>( TheFermionActionFactory::Instance().createObject(p.fermact.id, fermacttop, p.fermact.path) );
  }
  QDPIO::cout << "Creating State" << std::endl;
  state = S_asymm->createState(u);
  M_asymm =dynamic_cast<LinOpAsymm_T *>(S_asymm->linOp(state));


  // Do tests here

  T rhs= zero;
  gaussian(rhs,rb[1]);
  int success = 0;

  for(int i=0 ; i < p.iters; i++) {
    QDPIO::cout << "Doing test " << i+1 << " of " << p.iters << std::endl;
    // Zero the initial guesses
    Handle<MdagMMultiSystemSolver<T>> solver = dynamic_cast<MdagMMultiSystemSolver<T>*>(S_asymm->mInvMdagM(state,p.inv_param));

    auto n_shift = p.shifts.size();
    multi1d<T> solns(n_shift);
    for(int shift=0; shift < n_shift; ++shift) {
      solns[shift] = zero; // Zero all lattice
      //(solns[shift])[rb[1]]=zero;
    }

    (*solver)(solns,p.shifts,rhs);

    for(int shift = 0; shift < n_shift; ++shift) {
      T r;
      r[rb[1]] = zero;
      T tmp;
      tmp[rb[1]] = zero;

      (*M_asymm)(tmp,solns[shift],PLUS);
      (*M_asymm)(r, tmp, MINUS);

      r[rb[1]] += p.shifts[shift]*solns[shift];

      // -residudum
      r[rb[1]] -= rhs;

      Double resid_rel = sqrt(norm2(r,rb[1])/norm2(rhs,rb[1]));
      QDPIO::cout << "shift="<<shift << " || r || / || b ||=" << resid_rel << "   ";
      if ( toDouble(resid_rel) >= 1.0e-8 ) { 
        QDPIO::cout << "FAILED " << std::endl;
        success--; 
      }
      else { 
        QDPIO::cout << "PASSED " << std::endl;
      }
    }
  }

  // Done with tests
  Chroma::finalize();
  exit(success);
}



