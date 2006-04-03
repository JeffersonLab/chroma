// $Id: t_conslinop.cc,v 3.0 2006-04-03 04:59:14 edwards Exp $

#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,2};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml("t_conslinop.xml");
  push(xml,"t_conslinop");

  push(xml,"lattis");
  write(xml,"Nd",  Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);

  //! Example of calling a plaquette routine
  /*! NOTE: the STL is *not* used to hold gauge fields */
  multi1d<LatticeColorMatrix> u(Nd);

  QDPIO::cout << "Start gaussian\n";
  for(int m=0; m < u.size(); ++m)
    gaussian(u[m]);

  LatticeFermion psi, chi;
  gaussian(psi);

  {
    // Create a fermion BC. Note, the handle is on an ABSTRACT type
    Handle< FermState<LatticeFermion,
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> > > state(new PeriodicFermState<LatticeFermion,
					    multi1d<LatticeColorMatrix>,
					    multi1d<LatticeColorMatrix> >(u));

    WilsonDslash D(state);

    QDPIO::cout << "before dslash call" << endl;
    D.apply(chi, psi, PLUS, 0); 
    D.apply(chi, psi, PLUS, 1); 
    QDPIO::cout << "after dslash call" << endl;

    QDPIO::cout << "before wilson construct" << endl;
    Real Mass = 0.1;
    UnprecWilsonLinOp M(state,Mass);
    QDPIO::cout << "after wilson construct" << endl;
    M(chi, psi, PLUS); 
    QDPIO::cout << "after wilson call" << endl;
  }

  {  
    // Create your creator. Note, the handle is on an ABSTRACT type
    Handle< CreateFermState<LatticeFermion,
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> > >  cfs(new CreatePeriodicFermState<LatticeFermion,
					   multi1d<LatticeColorMatrix>,
					   multi1d<LatticeColorMatrix> >());

    Real Mass = 0.1;
    UnprecWilsonFermAct S(cfs,Mass);

    Handle< FermState<LatticeFermion,
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> > > state(S.createState(u));

    Handle< LinearOperator<LatticeFermion> > A(S.linOp(state));

    LatticeFermion   tmp;
    (*A)(tmp, psi, PLUS);
    DComplex np = innerProduct(psi,tmp);
    (*A)(tmp, psi, MINUS);
    DComplex nm = innerProduct(psi,tmp);

    push(xml,"norm_check");
    write(xml,"np", np);
    write(xml,"nm", nm);
    pop(xml);
  }

  pop(xml);

  // Time to bolt
  Chroma::finalize();

  exit(0);
}
