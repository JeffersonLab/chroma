// $Id: t_precact_sse.cc,v 1.1 2004-08-01 20:27:27 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_sse_w.h"

#include "qdp_util.h"

using namespace QDP;


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_precact.xml");
  push(xml, "t_precact");

  // Init the gauge field
  multi1d<LatticeColorMatrix> u(Nd);
//  HotSt(u);
  u = 1.0;

  InvType invType = CG_INVERTER;
  Real RsdCG = 1.0e-5;
  int MaxCG = 1000;
  int n_count = 0;

  // Create the BC objects
  const int bnd[] = {1,1,1,-1};
  multi1d<int> boundary(Nd);
  boundary = bnd;

  {
    // Create a fermion BC. Note, the handle is on an ABSTRACT type
    Handle< FermBC<multi1d<LatticeFermion> > >  fbc(new SimpleFermBC<multi1d<LatticeFermion> >(boundary));
    
    // The standard DWF fermact
    Real WilsonMass = 1.0;
    int N5 = 8;
    Real m_q = 0.3;
    EvenOddPrecDWFermActArray S_pdwf(fbc,WilsonMass,m_q,N5);
    Handle<const ConnectState> state(S_pdwf.createState(u));
    Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > A_pdwf(S_pdwf.linOp(state));
  
    SSEEvenOddPrecDWFermActArray S_sdwf(fbc,WilsonMass,m_q,N5);
    Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > A_sdwf(S_sdwf.linOp(state));
  
    // try the qprop
    multi1d<LatticeFermion>  psi5a(N5), psi5b(N5), chi5(N5), tmp1(N5);
    for(int m=0; m < N5; ++m)
    {
      gaussian(chi5[m]);
//      random(psi5a[m]);
      psi5a[m] = zero;
    }
    psi5b = psi5a;

    QDPIO::cout << "Prec inverter" << endl;
    S_pdwf.qpropT(psi5a, state, chi5, invType, RsdCG, MaxCG, n_count);
    QDPIO::cout << "SSE prec inverter" << endl;
    S_sdwf.opt_qpropT(psi5b, state, chi5, invType, RsdCG, MaxCG, n_count);
    
    for(int m=0; m < N5; ++m)
      tmp1[m] = psi5a[m] - psi5b[m];

    QDPIO::cout << "Test eo-prec and SSE opt eo-prec DWF qpropT" << endl
		<< "|pDWF|^2 = " << norm2(psi5a) << endl
		<< "|sDWF|^2 = " << norm2(psi5b) << endl
		<< "|pDWF - sDWF|^2 = " << norm2(tmp1) << endl;

#if 0
    // try the qpropT
    LatticeFermion chi, psia, psib;
    gaussian(chi);
    random(psia);
    psib = psia;

    QDPIO::cout << "Prec inverter" << endl;
    S_pdwf.qprop(psia, state, chi, invType, RsdCG, MaxCG, n_count);
    QDPIO::cout << "SSE prec inverter" << endl;
    S_sdwf.qprop(psib, state, chi, invType, RsdCG, MaxCG, n_count);
    
    QDPIO::cout << "Test eo-prec and SSE opt eo-prec DWF qprop" << endl
		<< "|pDWF|^2 = " << norm2(psia) << endl
		<< "|sDWF|^2 = " << norm2(psib) << endl
		<< "|pDWF - sDWF|^2 = " << norm2(psia-psib) << endl;
#endif
  }

  QDPIO::cout << "\n\n\n" << endl;

  pop(xml);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
