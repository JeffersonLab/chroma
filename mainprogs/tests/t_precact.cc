// $Id: t_precact.cc,v 1.8 2005-01-02 05:21:11 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

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
  HotSt(u);

  InvertParam_t  invParam;
  invParam.invType = CG_INVERTER;
  invParam.RsdCG = 1.0e-7;
  invParam.MaxCG = 1000;
  int n_count;

  // Create the BC objects
  const int bnd[] = {1,1,1,-1};
  multi1d<int> boundary(Nd);
  boundary = bnd;

  // Check a preconditioned and unpreconditioned fermact give the same results
  {
    // Create a fermion BC. Note, the handle is on an ABSTRACT type
    Handle< FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(boundary));

    // The Wilson fermact
    Real WilsonMass = -1;
    UnprecWilsonFermAct S_uwil(fbc, WilsonMass);
    Handle<const ConnectState> state(S_uwil.createState(u));
    Handle<const LinearOperator<LatticeFermion> > A_uwil(S_uwil.linOp(state));
  
    EvenOddPrecWilsonFermAct S_pwil(fbc, WilsonMass);
    Handle<const EvenOddPrecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> > > A_pwil(S_pwil.linOp(state));
  
    LatticeFermion  psi, chi, tmp1, tmp2;
    random(psi);
    QDPIO::cout << "Test unprec and eo-prec Wilson: sign=PLUS" << endl;
    (*A_uwil)(tmp1, psi, PLUS);
    (*A_pwil).unprecLinOp(tmp2, psi, PLUS);
    QDPIO::cout << "|Wil|^2 = " << norm2(tmp1) << endl
		<< "|pWil|^2 = " << norm2(tmp2) << endl
		<< "|Wil - pWil|^2 = " << norm2(tmp2-tmp1) << endl;

    QDPIO::cout << "Test unprec and eo-prec Wilson: sign=MINUS" << endl;
    (*A_uwil)(tmp1, psi, MINUS);
    A_pwil->unprecLinOp(tmp2, psi, MINUS);
    QDPIO::cout << "|Wil|^2 = " << norm2(tmp1) << endl
		<< "|pWil|^2 = " << norm2(tmp2) << endl
		<< "|Wil - pWil|^2 = " << norm2(tmp2-tmp1) << endl;

    // try the qprop
    gaussian(chi);
    random(tmp1);
    tmp2 = tmp1;
    QDPIO::cout << "Unprec Wilson inverter" << endl;
    {
      Handle<const SystemSolver<LatticeFermion> > qprop(S_uwil.qprop(state, invParam));
      (*qprop)(tmp1, chi);
    }
    QDPIO::cout << "Prec Wilson inverter" << endl;
    {
      Handle<const SystemSolver<LatticeFermion> > qprop(S_pwil.qprop(state, invParam));
      (*qprop)(tmp2, chi);
    }
    
    QDPIO::cout << "Test unprec and eo-prec Wilson inverter" << endl
		<< "|Wil|^2 = " << norm2(tmp1) << endl
		<< "|pWil|^2 = " << norm2(tmp2) << endl
		<< "|Wil - pWil|^2 = " << norm2(tmp2-tmp1) << endl;
  }

  QDPIO::cout << "\n\n\n" << endl;

  {
    // Create a fermion BC. Note, the handle is on an ABSTRACT type
    Handle< FermBC<multi1d<LatticeFermion> > >  fbc(new SimpleFermBC<multi1d<LatticeFermion> >(boundary));
    
    // The standard DWF fermact
    Real WilsonMass = 1.5;
    int N5 = 8;
    Real m_q = 0.1;
    UnprecDWFermActArray S_udwf(fbc,WilsonMass,m_q,N5);
    Handle<const ConnectState> state(S_udwf.createState(u));
    Handle<const LinearOperator< multi1d<LatticeFermion> > > A_udwf(S_udwf.linOp(state));
  
    EvenOddPrecDWFermActArray S_pdwf(fbc,WilsonMass,m_q,N5);
    Handle<const EvenOddPrecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> > > A_pdwf(S_pdwf.linOp(state));
  
    multi1d<LatticeFermion>  psi(N5), chi(N5), tmp1(N5), tmp2(N5);
    for(int m=0; m < N5; ++m)
      random(psi[m]);

    QDPIO::cout << "Test unprec and eo-prec DWF: sign=PLUS" << endl;
    (*A_udwf)(tmp1, psi, PLUS);
    A_pdwf->unprecLinOp(tmp2, psi, PLUS);
    QDPIO::cout << "|DWF|^2 = " << norm2(tmp1) << endl
		<< "|pDWF|^2 = " << norm2(tmp2) << endl;
    for(int m=0; m < N5; ++m)
      chi[m] = tmp2[m] - tmp1[m];
    QDPIO::cout << "|DWF - pDWF|^2 = " << norm2(chi) << endl;

    QDPIO::cout << "Test unprec and eo-prec DWF: sign=MINUS" << endl;
    (*A_udwf)(tmp1, psi, MINUS);
    A_pdwf->unprecLinOp(tmp2, psi, MINUS);
    QDPIO::cout << "|DWF|^2 = " << norm2(tmp1) << endl
		<< "|pDWF|^2 = " << norm2(tmp2) << endl;
    for(int m=0; m < N5; ++m)
      chi[m] = tmp2[m] - tmp1[m];
    QDPIO::cout << "|DWF - pDWF|^2 = " << norm2(chi) << endl;

    // try the qprop
    LatticeFermion chi5, psi5a, psi5b;
    gaussian(chi5);
    random(psi5a);
    psi5b = psi5a;

    QDPIO::cout << "Unprec inverter" << endl;
    {
      Handle<const SystemSolver<LatticeFermion> > qprop(S_udwf.qprop(state, invParam));
      (*qprop)(psi5a, chi5);
    }
    QDPIO::cout << "Prec inverter" << endl;
    {
      Handle<const SystemSolver<LatticeFermion> > qprop(S_pdwf.qprop(state, invParam));
      (*qprop)(psi5b, chi5);
    }
    
    QDPIO::cout << "Test unprec and eo-prec DWF inverter" << endl
		<< "|DWF|^2 = " << norm2(psi5a) << endl
		<< "|pDWF|^2 = " << norm2(psi5b) << endl
		<< "|DWF - pDWF|^2 = " << norm2(psi5a-psi5b) << endl;
  }

  QDPIO::cout << "\n\n\n" << endl;

  pop(xml);

  // Time to bolt
  QDP_finalize();

  exit(0);
}
