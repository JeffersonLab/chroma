// $Id: t_precact.cc,v 1.2 2003-11-25 03:03:56 edwards Exp $

#include <iostream>
#include <cstdio>

#include "chroma.h"

#include "qdp_util.h"


//---------------- HACK ----------------------------
// WARNING - Inefficient; improve later - move into QDP
#if 1
namespace QDP {

template<class T>
inline typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t
norm2(const multi1d< OLattice<T> >& s1)
{
  typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t  d;

  d = norm2(s1[0]);
  for(int n=1; n < s1.size(); ++n)
    d += norm2(s1[n]);

  return d;
}
}
#endif
//---------------------------------------------------


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

  InvType invType = CG_INVERTER;
  Real RsdCG = 1.0e-7;
  int MaxCG = 1000;
  int n_count;

  // Check a preconditioned and unpreconditioned fermact give the same results
  {
    // The Wilson fermact
    Real WilsonMass = -1;
    UnprecWilsonFermAct S_uwil(WilsonMass);
    const  LinearOperator<LatticeFermion>* A_uwil = S_uwil.linOp(u);
  
    EvenOddPrecWilsonFermAct S_pwil(WilsonMass);
    const  EvenOddPrecLinearOperator<LatticeFermion>* A_pwil = S_pwil.linOp(u);
  
    LatticeFermion  psi, chi, tmp1, tmp2;
    random(psi);
    QDPIO::cout << "Test unprec and eo-prec Wilson: sign=PLUS" << endl;
    (*A_uwil)(tmp1, psi, PLUS);
    A_pwil->unprecLinOp(tmp2, psi, PLUS);
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
    S_uwil.qprop(tmp1, u, chi, invType, RsdCG, MaxCG, n_count);
    QDPIO::cout << "Prec Wilson inverter" << endl;
    S_pwil.qprop(tmp2, u, chi, invType, RsdCG, MaxCG, n_count);
    
    QDPIO::cout << "Test unprec and eo-prec Wilson inverter" << endl
		<< "|Wil|^2 = " << norm2(tmp1) << endl
		<< "|pWil|^2 = " << norm2(tmp2) << endl
		<< "|Wil - pWil|^2 = " << norm2(tmp2-tmp1) << endl;
  }

  QDPIO::cout << "\n\n\n" << endl;

  {
    // The standard DWF fermact
    Real WilsonMass = 1.5;
    int N5 = 8;
    Real m_q = 0.1;
    UnprecDWFermActArray S_udwf(WilsonMass,m_q,N5);
    const  LinearOperator< multi1d<LatticeFermion> >* A_udwf = S_udwf.linOp(u);
  
    EvenOddPrecDWFermActArray S_pdwf(WilsonMass,m_q,N5);
    const  EvenOddPrecLinearOperator< multi1d<LatticeFermion> >* A_pdwf = S_pdwf.linOp(u);
  
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
    S_udwf.qprop(psi5a, u, chi5, invType, RsdCG, MaxCG, n_count);
    QDPIO::cout << "Prec inverter" << endl;
    S_pdwf.qprop(psi5b, u, chi5, invType, RsdCG, MaxCG, n_count);
    
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
