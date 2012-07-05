// $Id: t_precnef.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $

#include <iostream>
#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>

#include "chroma.h"

using namespace Chroma;

typedef multi1d<LatticeFermion>  MLF;
typedef multi1d<LatticeColorMatrix>  LCM;

typedef  void (EvenOddPrecLinearOperator< MLF, LCM >::* EO_mem)(MLF&, const MLF&, enum PlusMinus) const;


double time_func(const EvenOddPrecLinearOperator< MLF, LCM > *p, EO_mem A,
		 MLF& chi, const MLF& psi,
		 enum PlusMinus isign)
{
  clock_t myt1, myt2;
  double  mydt;
  int iter = 1;

  for(iter=1; ; iter <<= 1)
  {
    QDPIO::cout << "Applying D " << iter << " times" << endl;

    myt1=clock();
    for(int i=iter; i-- > 0; )
      (p->*A)(chi, psi, isign);
    myt2=clock();

    mydt=double(myt2-myt1)/double(CLOCKS_PER_SEC);
    QDPInternal::globalSum(mydt);
    mydt /= Layout::numNodes();

    if (mydt > 1)
      break;
  }

  myt1=clock();
  for(int i=iter; i-- > 0; )
    (p->*A)(chi, psi, isign);
  myt2=clock();

  mydt = (double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
  mydt *= 1.0e6/((double)(iter*(Layout::sitesOnNode()/2)));
  QDPInternal::globalSum(mydt);
  mydt /= Layout::numNodes();
  return mydt;
}


int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  //! Create a linear operator
  QDPIO::cout << "Constructing DWDslash" << endl;

  // Create a FermBC with only periodic BC. Note the handle is on an abstract type.
  Handle<FermBC<MLF> >  fbc_a(new PeriodicFermBC<MLF>);

  // DWDslash class can be optimised
  int N5 = 8;
  Real WilsonMass = 1.5;
  Real m_q = 0.01;

  Real b5=Real(1);
  Real c5=Real(0);
  EvenOddPrecNEFFermActArrayParams p;
  p.OverMass = WilsonMass;
  p.b5 = b5;
  p.c5 = c5;
  p.Mass = m_q;
  p.N5=N5;

  EvenOddPrecNEFFermActArray S_pdwf(fbc_a,p);

  Handle<const ConnectState> state(S_pdwf.createState(u));
  const EvenOddPrecLinearOperator< MLF, LCM >* D_pdwf = S_pdwf.linOp(state); 

  QDPIO::cout << "Done" << endl;

  MLF psi(S_pdwf.size()), chi(S_pdwf.size());
  for(int n=0; n < S_pdwf.size(); ++n)
    random(psi[n]);
  chi = zero;

  for(int isign = 1; isign >= -1; isign -= 2) 
  {
    QDPIO::cout << "Applying D" << endl;
    QDPIO::cout << " isign = " << isign << endl;
      
    PlusMinus is = (isign == 1 ? PLUS : MINUS);
    clock_t myt1;
    clock_t myt2;
    double mydt;

    int Ndiag  = 6*N5*Nc*Ns; // This is my count with the blas / chiral proj ops
    int NdiagInv = (10*N5-8)*Nc*Ns;
    int Neo    = 6*N5*Nc*Ns + N5*(1320+24);
    int Nflops = 2*Ndiag + 2*Neo + N5*24;

    // even-even-inv piece
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<MLF,LCM>::evenEvenInvLinOp, chi, psi, is);
    QDPIO::cout << "EvenEvenInv: The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(NdiagInv)/mydt) << ") Mflops " << endl;

    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<MLF,LCM>::evenEvenLinOp, chi, psi, is);
    QDPIO::cout << "EvenEven: The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Ndiag)/mydt) << ") Mflops " << endl;
      
    // odd-odd piece
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<MLF,LCM>::oddOddLinOp, chi, psi, is);
    QDPIO::cout << "OddOdd: The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Ndiag)/mydt) << ") Mflops " << endl;
      
    // even-odd
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<MLF,LCM>::evenOddLinOp, chi, psi, is);
    QDPIO::cout << "EvenOdd: The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Neo)/mydt) << ") Mflops " << endl;
    // odd-even
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<MLF,LCM>::oddEvenLinOp, chi, psi, (isign == 1 ? PLUS : MINUS));
    QDPIO::cout << "Odd-Even: The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Neo)/mydt) << ") Mflops " << endl;

    // Total thing
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<MLF,LCM>::operator(), chi, psi, is);
    QDPIO::cout << "Total: The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Nflops)/mydt) << ") Mflops " << endl;
  }
  
  delete D_pdwf;

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
