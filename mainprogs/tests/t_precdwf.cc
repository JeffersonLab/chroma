// $Id: t_precdwf.cc,v 1.3 2003-12-04 21:52:15 edwards Exp $

#include <iostream>
#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>

#include "chroma.h"

using namespace QDP;

typedef  void (EvenOddPrecLinearOperator< multi1d<LatticeFermion> >::* EO_mem)(multi1d<LatticeFermion>&, const multi1d<LatticeFermion>&, enum PlusMinus) const;


double time_func(const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > *p, EO_mem A,
		 multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi,
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
  const int foo[] = {8,4,4,4};
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
  multi1d<LatticeFermion> psi(N5), chi(N5);
  for(int n=0; n < N5; ++n)
    random(psi[n]);
  chi = zero;

  //! Create a linear operator
  QDPIO::cout << "Constructing DWDslash" << endl;

  // DWDslash class can be optimised
  Real WilsonMass = 1.5;
  Real m_q = 0.01;
  EvenOddPrecDWFermActArray S_pdwf(WilsonMass,m_q,N5);
  const ConnectStateProxy state(S_pdwf.createState(u));
  const  EvenOddPrecLinearOperator< multi1d<LatticeFermion> >* D_pdwf = S_pdwf.linOp(state); 

  QDPIO::cout << "Done" << endl;

  for(int isign = 1; isign >= -1; isign -= 2) 
  {
    QDPIO::cout << "Applying D" << endl;
    QDPIO::cout << " isign = " << isign << endl;
      
    PlusMinus is = (isign == 1 ? PLUS : MINUS);
    clock_t myt1;
    clock_t myt2;
    double mydt;

    int Ndiag  = (N5-2)*(5*24) + 2*(8*24);
    int Neo    = N5*(1320+24);
    int Nflops = 2*Ndiag + 2*Neo + N5*24;

    // even-even-inv piece
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<multi1d<LatticeFermion> >::evenEvenInvLinOp, chi, psi, is);
    QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Ndiag)/mydt) << ") Mflops " << endl;
      
    // odd-odd piece
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<multi1d<LatticeFermion> >::oddOddLinOp, chi, psi, is);
    QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Ndiag)/mydt) << ") Mflops " << endl;
      
    // even-odd
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<multi1d<LatticeFermion> >::evenOddLinOp, chi, psi, is);
    QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Neo)/mydt) << ") Mflops " << endl;
    // odd-even
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<multi1d<LatticeFermion> >::oddEvenLinOp, chi, psi, (isign == 1 ? PLUS : MINUS));
    QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Neo)/mydt) << ") Mflops " << endl;

    // Total thing
    mydt = time_func(D_pdwf, &EvenOddPrecLinearOperator<multi1d<LatticeFermion> >::operator(), chi, psi, is);
    QDPIO::cout << "The time per lattice point is "<< mydt << " micro sec" 
		<< " (" <<  ((double)(Nflops)/mydt) << ") Mflops " << endl;
  }
  
  delete D_pdwf;

  // Time to bolt
  QDP_finalize();

  exit(0);
}
