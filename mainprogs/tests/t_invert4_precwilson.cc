// $Id: t_invert4_precwilson.cc,v 3.1 2009-09-25 12:41:23 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"

#include "actions/ferm/invert/syssolver_linop_cg.h"
#include "actions/ferm/invert/syssolver_linop_cg.h"

#include "actions/ferm/invert/quda_solvers/syssolver_linop_cg_quda_wilson_single.h"

using namespace Chroma;


enum GaugeStartType { COLD_START=0, HOT_START };
struct Params_t { 
  multi1d<int> nrow;
  multi1d<int> boundary;
  GaugeStartType gauge_start_type;
};



void checkInverter(multi1d<LatticeColorMatrix>& u)
{
  LatticeFermion psi;
  LatticeFermion psi2;
  LatticeFermion chi;
  typedef LatticeFermion T;
  typedef multi1d<LatticeColorMatrix> Q;
  typedef multi1d<LatticeColorMatrix> P;


  multi1d<int> boundary(4);
  boundary[0] = 1; 
  boundary[1] = 1;
  boundary[2] = 1;
  boundary[3] = 1;
  Real mass = 0;

  // Set Up CG Solver
  SysSolverCGParams cg_p;
  cg_p.RsdCG = 1.0e-8;
  cg_p.MaxCG = 1000;
  cg_p.MinCG = 5;

  SysSolverCGQUDAWilsonParams quda_p;
  quda_p.WilsonParams.Mass = mass;
  quda_p.WilsonParams.anisoParam.anisoP = false;
  quda_p.WilsonParams.anisoParam.t_dir = 3;
  quda_p.WilsonParams.anisoParam.xi_0 = 1.0;
  quda_p.WilsonParams.anisoParam.nu = 1.0;
  quda_p.AntiPeriodicT = false;
  quda_p.MaxIter =cg_p.MaxCG;
  quda_p.RsdTarget = cg_p.RsdCG;
  quda_p.Delta = 1.0e-10;


  Handle<FermBC<T,P,Q> >  
    fbc(new SimpleFermBC<T,P,Q>(boundary));

   // Create a FermState Creator with boundaries
  Handle<CreateFermState<T,P,Q> > cfs( new CreateSimpleFermState<T,P,Q>(fbc));
  EvenOddPrecWilsonFermAct  S_w(cfs, mass);

  
  // Apply boundary to u
  Handle< FermState<T,P,Q> > connect_state(S_w.createState(u));

  Handle< LinearOperator<LatticeFermion> > D_op( S_w.linOp(connect_state) );

  // Get Initial Vector
  chi=zero;
  gaussian(chi,rb[1]);
  psi=zero;
  psi2=zero;

  QDPIO::cout << "chi_norm = " << sqrt(norm2(chi, rb[1])) << endl;

  
  LinOpSysSolverCG<LatticeFermion> CGSolver(D_op,cg_p);
  CGSolver(psi, chi);

  LinOpSysSolverCGQUDAWilson QUDA_CG_Solver(D_op, connect_state, quda_p);
  QUDA_CG_Solver(psi2,chi);

  LatticeFermion r;
  r[rb[1]] = psi2 - psi;
  Double chi_norm_diff = norm2(r, rb[1]);
  QDPIO::cout << " || chi2 - chi || = " << chi_norm_diff << endl;

}

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  multi1d<int> nrow(Nd);
  for(int i=0; i < Nd-1; i++) { 
    nrow[i] = 4;
  }
  nrow[Nd-1] = 16;
  Layout::setLattSize(nrow);
  Layout::create();

    struct Cfg_t config = { CFG_TYPE_WEAK_FIELD, "dummy" };
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, config);
  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Setup the lattice
 
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"t_invert");
  push(xml_out,"params");
  write(xml_out, "nrow", nrow);
  pop(xml_out); // Params

  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  checkInverter(u);

  pop(xml_out);
  xml_out.close();

  Chroma::finalize();
    
  exit(0);
}
