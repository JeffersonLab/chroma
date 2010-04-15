// $Id: t_invert4_precwilson.cc,v 3.4 2009-10-09 13:59:46 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"

#include "io/xml_group_reader.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include <string>

#include "actions/ferm/invert/minvcg2.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_cg_chrono_clover.h"
using namespace Chroma;
using namespace std;

struct AppParams { 
  multi1d<int> nrow;
  Cfg_t inputCfg;
  GroupXML_t fermact;
  MultiSysSolverCGChronoCloverParams invParam;
};



void checkInverter(const AppParams& p, multi1d<LatticeColorMatrix>& u)
{

  typedef LatticeFermion T;
  typedef multi1d<LatticeColorMatrix> Q;
  typedef multi1d<LatticeColorMatrix> P;

  typedef LatticeFermionF TF;
  typedef multi1d<LatticeColorMatrixF> QF;
  typedef multi1d<LatticeColorMatrixF> PF;

  std::istringstream is(p.fermact.xml);
  XMLReader fermact_xml(is);
  Handle< WilsonTypeFermAct<T,P,Q> > S_f ( TheWilsonTypeFermActFactory::Instance().createObject(p.fermact.id, fermact_xml, p.fermact.path) );

  
    //S_f( 
    //				 new EvenOddPrecCloverFermAct(CreateFermStateEnv::reader(fermact_xml, p.fermact.path), CloverFermActParams(fermact_xml, p.fermact.path)));
  
  
  Handle< FermState<T,P,Q> > connect_state( S_f->createState(u) );
  Handle< LinearOperator<LatticeFermion> > D_op( S_f->linOp(connect_state) );
  

  MdagMMultiSysSolverCGChronoClover mprec(D_op, connect_state, p.invParam);

  int n_shifts=12;
  multi1d<Real> shifts(n_shifts);
  shifts[0]  =  1.51480322230288e-05  ;
  shifts[1]  =  0.000165259114300968  ;
  shifts[2]  =  0.000661184012887753  ;
  shifts[3]  =  0.00214955953261573  ;
  shifts[4]  =  0.00657190811374048 ;
  shifts[5]  =  0.0197043312297316   ;
  shifts[6]  =  0.0587640674448055  ;
  shifts[7]  =  0.175534749714768  ;
  shifts[8]  =  0.529999019144126  ;
  shifts[9]  =  1.6565878795396  ;
  shifts[10] = 5.78532483906374   ;
  shifts[11] = 30.6835527589452  ;

  multi1d<Real> RsdCG(shifts.size());
  if (p.invParam.RsdTarget.size() == 1) {
    RsdCG = p.invParam.RsdTarget[0];
  }
  else if (p.invParam.RsdTarget.size() == RsdCG.size()) {
    
    RsdCG = p.invParam.RsdTarget;
  }
  else {
    
    QDPIO::cerr << "MdagMMultiSysSolverCGChronoClover: shifts incompatible" << endl;
    QDP_abort(1);
  }


  
  LatticeFermion chi; 
  gaussian(chi, D_op->subset());

  // Initialized psi
  multi1d<LatticeFermion> psi(n_shifts);
  for(int i=0; i < n_shifts; i++){ 
    psi[i] = zero;
  }

  SystemSolverResults_t res;
  StopWatch swatch;
  swatch.reset();
  swatch.start();
  MInvCG2(*D_op, chi, psi, shifts, RsdCG, p.invParam.MaxIter, res.n_count);
  swatch.stop();

  // Check solutions
  for(int i=0; i < n_shifts; i++) { 
    LatticeFermion r;
    r[D_op->subset()]= chi;
    LatticeFermion tmp1,tmp2;
    (*D_op)(tmp1, psi[i], PLUS);
    (*D_op)(tmp2, tmp1, MINUS);
    tmp2[D_op->subset()] += shifts[i]*psi[i];
    r[D_op->subset()] -= tmp2;
    Double norm = norm2(r, D_op->subset())/norm2(chi, D_op->subset());

    QDPIO::cout << "soln "<< i <<" : RsdCG=" << RsdCG[i] << "  r=" << sqrt(norm) << endl;
  }
  QDPIO::cout << "MinvCG: n_count = " << res.n_count << "  Time=" << swatch.getTimeInSeconds() << endl; 

  swatch.reset();
  swatch.start();

  SystemSolverResults_t res2 = mprec(psi, shifts, chi);
    
  swatch.stop();
  for(int i=0; i < n_shifts; i++) { 
    LatticeFermion r;
    r[D_op->subset()]= chi;
    LatticeFermion tmp1,tmp2;
    (*D_op)(tmp1, psi[i], PLUS);
    (*D_op)(tmp2, tmp1, MINUS);
    tmp2[D_op->subset()] += shifts[i]*psi[i];
    r[D_op->subset()] -= tmp2;
    Double norm = norm2(r, D_op->subset())/norm2(chi, D_op->subset());

    QDPIO::cout << "soln "<< i <<" : RsdCG=" << RsdCG[i] << "  r=" << sqrt(norm) << endl;
  }
  QDPIO::cout << "Chrono: n_count = " << res2.n_count << "  Time=" << swatch.getTimeInSeconds() << endl; 

}


void read(XMLReader& r, const std::string path, AppParams& p) 
{
  read(r, "nrow", p.nrow);
  read(r, "Cfg", p.inputCfg);
  p.fermact = readXMLGroup(r, "FermionAction", "FermAct");
  read(r, "InvertParam", p.invParam);

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
