// $Id: t_precact_sse.cc,v 1.3 2004-09-09 15:52:52 edwards Exp $

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
  const int foo[] = {4,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_precact.xml");
  push(xml, "t_precact");

  // Init the gauge field
  multi1d<LatticeColorMatrix> u(Nd);
//  HotSt(u);
//  u = 1.0;
//  XMLReader gauge_xml;
//  readSzin(gauge_xml, u, "test_purgaug.cfg1");
  for(int m=0; m < Nd; ++m)
  {
    ColorMatrix t;
    t = 0;
    for(int i=0; i < Nc; ++i)
      for(int j=0; j < Nc; ++j)
	pokeColor(t,cmplx(Real((m+1)*((i+j)*0.02)),Real(-(m+1)*((i+j)*0.00))),i,j);

#if 1
    for(int site=0; site < Layout::vol(); ++site)
      pokeSite(u[m],ColorMatrix(site+t),Layout::siteCoords(0,site));
#else
    u[m] = t;
#endif
//    reunit(u[m]);
  }

  InvertParam_t  invParam;
  invParam.invType = CG_INVERTER;
  invParam.RsdCG = 1.0e-12;
  invParam.MaxCG = 3000;
  int n_count = 0;

  // Create the BC objects
  const int bnd[] = {1,1,1,1};
  multi1d<int> boundary(Nd);
  boundary = bnd;

  {
    // Create a fermion BC. Note, the handle is on an ABSTRACT type
    Handle< FermBC<multi1d<LatticeFermion> > >  fbc(new SimpleFermBC<multi1d<LatticeFermion> >(boundary));
    
    // The standard DWF fermact
    Real WilsonMass = 1.1;
    int N5 = 8;
    Real m_q = 0.3;
    UnprecDWFermActArray S_pdwf(fbc,WilsonMass,m_q,N5);
    Handle<const ConnectState> state(S_pdwf.createState(u));
  
    SSEEvenOddPrecDWFermActArray S_sdwf(fbc,WilsonMass,m_q,N5);
  
    multi1d<int> coord(Nd);
    coord = 0;

    // try the qprop
    multi1d<LatticeFermion>  psi5a(N5), psi5b(N5), chi5(N5), tmp1(N5);
    for(int m=0; m < N5; ++m)
    {
      gaussian(chi5[m]);
      random(psi5a[m]);
//      srcfil(chi5[m], coord, 0, 0);
//      psi5a[m] = zero;
    }
    psi5b = psi5a;

    QDPIO::cout << "UnPrec inverter" << endl;
    S_pdwf.qpropT(psi5a, state, chi5, invParam, n_count);
    QDPIO::cout << "SSE prec inverter" << endl;
    S_sdwf.opt_qpropT(psi5b, state, chi5, invParam, n_count);
    
    for(int m=0; m < N5; ++m)
      tmp1[m] = psi5a[m] - psi5b[m];

    QDPIO::cout << "Test eo-prec and SSE opt eo-prec DWF qpropT" << endl
		<< "|pDWF|^2 = " << norm2(psi5a) << endl
		<< "|sDWF|^2 = " << norm2(psi5b) << endl
		<< "|pDWF - sDWF|^2 = " << norm2(tmp1) << endl;

    
    push(xml,"Unprec");
    write(xml,"psi5a",psi5a);
    pop(xml);

    push(xml,"Prec");
    write(xml,"psi5b",psi5b);
    pop(xml);


#if 0
    // try the qpropT
    LatticeFermion chi, psia, psib;
    gaussian(chi);
    random(psia);
    psib = psia;

    QDPIO::cout << "UnPrec inverter" << endl;
    S_pdwf.qprop(psia, state, chi, invParam, n_count);
    QDPIO::cout << "SSE prec inverter" << endl;
    S_sdwf.qprop(psib, state, chi, invParam, n_count);
    
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
