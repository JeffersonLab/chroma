// $Id: t_invert4_precwilson.cc,v 3.4 2009-10-09 13:59:46 bjoo Exp $

#include "chroma.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>


#include "io/xml_group_reader.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include <string>

#include "actions/ferm/invert/minvcg2.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"
#include "actions/ferm/invert/mg_proto/mgproto_solver_params.h"

#include "actions/ferm/invert/mg_proto/syssolver_linop_clover_mg_proto.h"
#include "actions/ferm/invert/mg_proto/mg_proto_helpers.h"

#include "lattice/fine_qdpxx/wilson_clover_linear_operator.h"

using namespace Chroma;

struct AppParams { 
  multi1d<int> nrow;
  Cfg_t inputCfg;
  GroupXML_t fermact;
  MGProtoSolverParams mgproto_params;
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
  Handle< Chroma::LinearOperator<LatticeFermion> > D_op( S_f->linOp(connect_state) );
  

  XMLBufferWriter xml_buf;
  write(xml_buf, "MGProtoParams", p.mgproto_params);
  QDPIO::cout << xml_buf.str() << std::endl;


  // Check the FineM against the linear operator.
  QDPIO::cout << "Testing MGProto Fine Linear Operator: " << std::endl;
  {
	  LatticeFermion psi; LatticeFermion chi1; LatticeFermion chi2;
	  gaussian(psi);
	  chi1=zero;
	  chi2=zero;


	  QDPIO::cout << "Applying Chroma UnprecCloverLinop" <<std::endl;
	  (*D_op)(chi1,psi,PLUS);


	  std::shared_ptr<const MG::QDPWilsonCloverLinearOperator> D_op_mg = MGProtoHelpers::createFineLinOp(p.mgproto_params,
			  connect_state->getLinks());

	  QDPIO::cout << "Applying MG QDPWilsonCloverLinearOperator " << std::endl;
	  (*D_op_mg)(chi2,psi,MG::LINOP_OP);

	  LatticeFermion diff=chi2-chi1;
	  Double n2_diff = sqrt(norm2(diff));
	  QDPIO::cout << "Diff Norm = " << n2_diff << std::endl;
	  QDPIO::cout << "Diff Norm/norm2psi = " << n2_diff /sqrt(norm2(psi)) << std::endl;
	  QDPIO::cout << "Diff Norm / d.o.f = " << n2_diff / toDouble(Layout::vol()*Nc*Ns) << std::endl;
  }

  // Try and create the solver object:

  {
	  LinOpSysSolverMGProtoClover  the_solver(D_op, connect_state, p.mgproto_params);
	  LatticeFermion psi, chi;
	  gaussian(chi);
	  psi=zero;

	  the_solver(psi,chi);

	  LatticeFermion diff;
	  diff=zero;

	  (*D_op)(diff,psi, PLUS);   // Diff = Ax
	  diff -= chi;         // Diff = Ax - b = -( b - Ax )
	  Double n2diff = norm2(diff);
	  Double reln2diff = n2diff/norm2(chi);
	  QDPIO::cout << " || b - A x || = "  << sqrt(n2diff) << std::endl;
	  QDPIO::cout << " || b - A x || /  || b || = " << sqrt(reln2diff) << std::endl;

  }

  // Now do it again. This time one should fine the solver object
  {
 	  LinOpSysSolverMGProtoClover  the_solver(D_op, connect_state, p.mgproto_params);
 	  LatticeFermion psi, chi;
 	  gaussian(chi);
 	  psi=zero;

 	  the_solver(psi,chi);

 	  LatticeFermion diff;
 	  diff=zero;
 	  (*D_op)(diff,psi, PLUS);   // Diff = Ax
 	  diff -= chi;         // Diff = Ax - b = -( b - Ax )
 	  Double n2diff = norm2(diff);
 	  Double reln2diff = n2diff/norm2(chi);
 	  QDPIO::cout << " || b - A x || = "  << sqrt(n2diff) << std::endl;
 	  QDPIO::cout << " || b - A x || /  || b || = " << sqrt(reln2diff) << std::endl;

   }

  // Delete the MG Stuff
  MGProtoHelpers::deleteMGPreconditioner(p.mgproto_params.SubspaceId);


}


void read(XMLReader& r, const std::string path, AppParams& p) 
{
  read(r, "nrow", p.nrow);
  read(r, "Cfg", p.inputCfg);
  p.fermact = readXMLGroup(r, "FermionAction", "FermAct");
  read(r, "MGProtoParams", p.mgproto_params);
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
  QDPIO::cout << "Linkage = " << linkageHack() << std::endl;


  AppParams params;

  XMLReader xml_in;
  try
  {
    xml_in.open(Chroma::getXMLInputFileName());
    read(xml_in, "/Param", params);
  }
  catch(const std::string& e) 
  {
    QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
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
