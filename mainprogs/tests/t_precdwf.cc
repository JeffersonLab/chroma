// $Id: t_precdwf.cc,v 1.10 2005-03-01 19:12:53 edwards Exp $

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
    Internal::broadcast(mydt);

    if (mydt > 1)
      break;
  }

  myt1=clock();
  for(int i=iter; i-- > 0; )
    (p->*A)(chi, psi, isign);
  myt2=clock();

  mydt = (double)(myt2-myt1)/((double)(CLOCKS_PER_SEC));
  mydt *= 1.0e6/((double)(iter*(Layout::sitesOnNode()/2)));
  return mydt;
}


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {8,8,8,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_lwldslash.xml");

  //! Test out dslash
  multi1d<LatticeColorMatrix> u(Nd);
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  int N5 = 8;
  MLF psi(N5), chi(N5);
  for(int n=0; n < N5; ++n)
    random(psi[n]);
  chi = zero;

  //! Create a linear operator
  QDPIO::cout << "Constructing DWDslash" << endl;

  // Create a FermBC with only periodic BC. Note the handle is on an abstract type.
  Handle<FermBC<MLF> >  fbc_a(new PeriodicFermBC<MLF>);

  // DWDslash class can be optimised
  Real WilsonMass = 1.5;
  Real m_q = 0.01;
  EvenOddPrecDWFermActArray S_pdwf(fbc_a,WilsonMass,m_q,N5);
  Handle<const ConnectState> state(S_pdwf.createState(u));
  const EvenOddPrecLinearOperator< MLF, LCM >* D_pdwf = S_pdwf.linOp(state); 

  QDPIO::cout << "Done" << endl;

  for(int isign = 1; isign >= -1; isign -= 2) 
  {
    QDPIO::cout << "Applying D" << endl;
    QDPIO::cout << " isign = " << isign << endl;
      
    PlusMinus is = (isign == 1 ? PLUS : MINUS);
    clock_t myt1;
    clock_t myt2;
    double mydt;

//    int Ndiag  = (N5-2)*(5*24) + 2*(8*24);
    int Ndiag  = N5*(4*24) + (N5-1)*(8*24) + 3*24;   // this is what I get counting flops in code
    int Neo    = N5*(1320+24);
    int Nflops = 2*Ndiag + 2*Neo + N5*24;

    // even-even-inv piece
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<MLF,LCM>::evenEvenInvLinOp, chi, psi, is);
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
  QDP_finalize();

  exit(0);
}
