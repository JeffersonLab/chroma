// $Id: t_precact_sse.cc,v 3.1 2006-07-03 15:26:11 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
// #include "actions/ferm/qprop/prec_dwf_qprop_array_sse_w.h"

#include "qdp_util.h"

using namespace Chroma;
using namespace Chroma;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_precact.xml");
  push(xml, "t_precact");

  // Init the gauge field
  multi1d<LatticeColorMatrix> u(Nd);
  HotSt(u);
//  u = 1.0;
//  XMLReader gauge_xml;
//  readSzin(gauge_xml, u, "test_purgaug.cfg1");

  SysSolverCGParams  invParam;
  invParam.RsdCG = 1.0e-6;
  invParam.MaxCG = 3000;
  int n_count = 0;

  GroupXML_t inv_param;
  {
    XMLBufferWriter xml_buf;
    write(xml_buf, "InvertParam", invParam);
    XMLReader xml_in(xml_buf);
    inv_param = readXMLGroup(xml_in, "/InvertParam", "invType");
  }

  // Create the BC objects
  const int bnd[] = {1,1,1,-1};
  multi1d<int> boundary(Nd);
  boundary = bnd;

  // Create a fermion BC. Note, the handle is on an ABSTRACT type
  Handle< FermBC<multi1d<LatticeFermion> > >  fbc(new SimpleFermBC<multi1d<LatticeFermion> >(boundary));
    
  // The standard DWF fermact
  Real WilsonMass = 1.1;
  int N5 = 8;
  Real m_q = 0.3;
  EvenOddPrecDWFermActArray S_pdwf(fbc,WilsonMass,m_q,N5);
  Handle<const ConnectState> state(S_pdwf.createState(u));
  
  {
    // Setup solvers
//    Handle< const SystemSolver< multi1d<LatticeFermion> > > qpropT(S_pdwf.qpropT(state,inv_param));
    Handle< const SystemSolver< multi1d<LatticeFermion> > > QDPqpropT(new PrecFermAct5DQprop<LatticeFermion, multi1d<LatticeColorMatrix> >(Handle< const EvenOddPrecLinearOperator<multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > >(S_pdwf.linOp(state)), inv_param));

    Handle< const SystemSolver< multi1d<LatticeFermion> > > SSEqpropT(new SSEDWFQpropT(state,WilsonMass,m_q,N5,inv_param));

    // Try the qpropT
    multi1d<LatticeFermion>  psi5a(N5), psi5b(N5), chi5(N5), tmp1(N5);
    for(int m=0; m < N5; ++m)
    {
      gaussian(chi5[m]);
      random(psi5a[m]);
    }
    psi5b = psi5a;

    QDPIO::cout << "Prec inverter: " << endl;
    QDPIO::cout << "  iterations = " << (*QDPqpropT)(psi5a, chi5) << endl;
    QDPIO::cout << "SSE prec inverter" << endl;
    QDPIO::cout << "  iterations = " << (*SSEqpropT)(psi5b, chi5) << endl;
    
    for(int m=0; m < N5; ++m)
      tmp1[m] = psi5a[m] - psi5b[m];

    QDPIO::cout << "Test eo-prec and SSE opt eo-prec DWF qpropT" << endl
		<< "|pDWF|^2 = " << norm2(psi5a) << endl
		<< "|sDWF|^2 = " << norm2(psi5b) << endl
		<< "|pDWF - sDWF|^2 = " << norm2(tmp1) << endl;
  }

  {
    // Setup solvers
    Handle< const SystemSolver< LatticeFermion > > qprop(S_pdwf.qprop(state,inv_param));

    // Try the qprop
    LatticeFermion chi, psia, psib;
    gaussian(chi);
    random(psia);
    psib = psia;

    QDPIO::cout << "Prec inverter" << endl;
    QDPIO::cout << "  iterations = " << (*qprop)(psia, chi) << endl;
    QDPIO::cout << "SSE prec inverter" << endl;
    QDPIO::cout << "  iterations = " << (*qprop)(psib, chi) << endl;
    
    QDPIO::cout << "Test eo-prec and SSE opt eo-prec DWF qprop" << endl
		<< "|pDWF|^2 = " << norm2(psia) << endl
		<< "|sDWF|^2 = " << norm2(psib) << endl
		<< "|pDWF - sDWF|^2 = " << norm2(psia-psib) << endl;
  }

  QDPIO::cout << "\n\n\nDone" << endl;

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
