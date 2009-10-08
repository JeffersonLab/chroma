// $Id: t_invert4_precwilson.cc,v 3.3 2009-10-08 00:52:23 bjoo Exp $

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
#include "io/xml_group_reader.h"
#include "actions/ferm/invert/quda_solvers/syssolver_quda_clover_params.h"
#include "actions/ferm/invert/quda_solvers/syssolver_linop_quda_clover.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include <string>
using namespace Chroma;
using namespace std;

struct AppParams { 
  multi1d<int> nrow;
  Cfg_t inputCfg;
  GroupXML_t fermact;
  GroupXML_t quda_solver;
};



void checkInverter(const AppParams& p, multi1d<LatticeColorMatrix>& u)
{

  typedef LatticeFermion T;
  typedef multi1d<LatticeColorMatrix> Q;
  typedef multi1d<LatticeColorMatrix> P;

  std::istringstream is(p.fermact.xml);
  XMLReader fermact_xml(is);
  Handle< EvenOddPrecCloverFermAct > S_f( 
					 new EvenOddPrecCloverFermAct(CreateFermStateEnv::reader(fermact_xml, p.fermact.path), CloverFermActParams(fermact_xml, p.fermact.path)));
  
  
  Handle< FermState<T,P,Q> > connect_state( S_f->createState(u) );

  Handle< LinearOperator<LatticeFermion> > D_op( S_f->linOp(connect_state) );

  std::istringstream is2(p.quda_solver.xml);
  XMLReader quda_xml(is2);

  Handle<  LinOpSystemSolver<LatticeFermion> > QUDASolver( new LinOpSysSolverQUDAClover(D_op, connect_state, SysSolverQUDACloverParams(quda_xml, p.quda_solver.path)));
							 
T chi, psi;

// Get Initial Vector
  chi=zero;
  gaussian(chi,rb[1]);
  psi=zero;
  (*QUDASolver)(psi,chi);

  T r = chi;
T tmp; 
(*D_op)(tmp, psi, PLUS);
r[rb[1]] -= tmp; 


QDPIO::cout << "||r||= " << sqrt(norm2(r, rb[1])) << endl;
QDPIO::cout << "||r||/||b||= " << sqrt(norm2(r, rb[1]))/sqrt(norm2(chi,rb[1])) << endl;

  
  
}


void read(XMLReader& r, const std::string path, AppParams& p) 
{
  read(r, "nrow", p.nrow);
  read(r, "Cfg", p.inputCfg);
  p.fermact = readXMLGroup(r, "FermionAction", "FermAct");
  p.quda_solver = readXMLGroup(r, "InvertParam", "invType");
}

bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= InlineAggregateEnv::registerAll();
  foo &= GaugeInitEnv::registerAll();

  return foo;
}

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);
  QDPIO::cout << "Linkage = " << linkageHack() << endl;


  AppParams params;

  XMLReader xml_in;
  try
  {
    xml_in.open(Chroma::getXMLInputFileName());
    read(xml_in, "/Param", params);
  }
  catch(const std::string& e) 
  {
    QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
    QDP_abort(1);
  }
  Layout::setLattSize(params.nrow);
  Layout::create();
  
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, params.inputCfg);
  unitarityCheck(u);

  // Setup the lattice
 
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"t_invert");
  push(xml_out,"params");
  write(xml_out, "nrow", params.nrow);
  pop(xml_out); // Params

  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  checkInverter(params, u);

  pop(xml_out);
  xml_out.close();

  Chroma::finalize();
    
  exit(0);
}
